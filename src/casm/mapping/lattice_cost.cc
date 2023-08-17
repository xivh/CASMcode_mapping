#include "casm/mapping/lattice_cost.hh"

#include "casm/crystallography/SymType.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace mapping {

/// \brief Returns the volume-normalized strain cost, calculated to be
/// invariant to which structure is the child/parent.
///
/// Given a lattice mapping:
/// \f[
///     F * L_1 * T_1 * N = L_2
/// \f]
///
/// Or equivalently,
/// \f[
///     L_1 * T_1 * N = F_reverse * L_2
/// \f]
///
/// This calculates the cost as:
///
/// Parent to child deformation cost:
/// \f[
///       \tilde{U} = \frac{1}{\det{U}^{1/3}} U
///       \tilde{B} = \tilde{U} - I
///       cost = (1./3.)*\mathrm{tr}(\tilde{B}^{2})
/// \f]
///
/// Child to parent deformation cost, using \f$V_{reverse}^{-1} = U\f$:
/// \f[
///       \tilde{V}_{reverse} = \frac{1}{\det{V_{reverse}}^{1/3}} V_{reverse}
///       \tilde{B}_{reverse} = \tilde{V_{reverse}} - I
///       cost = (1./3.)*\mathmr{Tr}{\tilde{B}_{reverse}^{2}}
/// \f]
///
/// Direction invariant cost:
/// \f[
///       strain_cost = (1./2.)*(
///           (1./3.)*\mathrm{tr}(\tilde{B}^{2}) +
///           (1./3.)*\mathrm{tr}(\tilde{B}_{reverse}^{2}))
/// \f]
///
/// In the above, the \f$\tilde\f$ indicates that the value is normalized to be
/// volume invariant.
///
/// \param deformation_gradient The deformation gradient, \f$F or F_reverse\f$.
///     The result is equivalent whether this is the parent to child
///     deformation or child to parent deformation.
///
double isotropic_strain_cost(Eigen::Matrix3d const &deformation_gradient) {
  Eigen::Matrix3d const &F = deformation_gradient;
  double vol_factor = std::pow(std::abs(F.determinant()), 1. / 3.);

  Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);
  Eigen::Matrix3d U_normalized = polar_decomposition(F / vol_factor);
  Eigen::Matrix3d V_reverse_normalized = U_normalized.inverse();

  // M.squaredNorm() = squared Frobenius norm of M = tr(M*M.transpose())
  return ((U_normalized - I).squaredNorm() +
          (V_reverse_normalized - I).squaredNorm()) /
         6.;
}

/// \brief Returns the symmetrized right stretch tensor
///
/// \param deformation_gradient The deformation gradient.
/// \param lattice1_point_group Reference lattice point group. Use the crystal
/// point group
///     of the prim if mapping structures.
///
/// \return \f$U_{symmetrized}\f$, where:
///   - \f$U_{symmetrized} = (1/N_{G_1}) * \sum_i (G_1(i) * U * G_1(i)^{-1})\f$,
///   - \f$U\f$: right stretch tensor of the deformation gradient
///   - \f$N_{G_1}\f$: Reference lattice point group size
///   - \f$G_1(i)\f$: Element \f$i\f$ of the reference lattice point group
///
Eigen::Matrix3d symmetrized_right_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group) {
  Eigen::Matrix3d const &F = deformation_gradient;
  Eigen::Matrix3d U = polar_decomposition(F);

  Eigen::Matrix3d U_aggregate = Eigen::Matrix3d::Zero();
  for (auto const &op : lattice1_point_group) {
    U_aggregate += op.matrix * U * op.matrix.inverse();
  }
  return U_aggregate / double(lattice1_point_group.size());
}

/// \brief Returns the symmetry-breaking component of the right stretch tensor
///
/// \param deformation_gradient The deformation gradient.
/// \param lattice1_point_group Reference lattice point group. Use the crystal
/// point group
///     of the prim if mapping structures.
///
/// \return \f$U - U_{symmetrized} + I\f$, where:
///   - \f$U\f$: Right stretch tensor of the deformation gradient
///   - \f$U_{symmetrized} = (1/N_{G_1}) * \sum_i (G_1(i) * U * G_1(i)^{-1})\f$,
///   - \f$N_{G_1}\f$: Parent point group size
///   - \f$G_1(i)\f$: Element \f$i\f$ of the reference lattice point group
///   - \f$I\f$: Identity matrix
///
Eigen::Matrix3d symmetry_breaking_right_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group) {
  Eigen::Matrix3d const &F = deformation_gradient;
  Eigen::Matrix3d U = polar_decomposition(F);

  Eigen::Matrix3d U_aggregate = Eigen::Matrix3d::Zero();
  for (auto const &op : lattice1_point_group) {
    U_aggregate += op.matrix * U * op.matrix.inverse();
  }

  Eigen::Matrix3d U_sym = U_aggregate / double(lattice1_point_group.size());

  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d U_sym_breaking = U - U_sym + I;
  return U_sym_breaking;
}

/// \brief Returns the symmetry-breaking strain cost, calculated to be
/// invariant to which structure is the child/parent.
///
/// The symmetry-breaking strain cost is calculated:
///
/// \f[
///     \frac{1}{2}*\left(
///         \frac{1}{3}*\mathrm{tr}\left( B_{sym-break}^2 \right) +
///         \frac{1}{3}*\mathrm{tr}\left( B_{reverse-sym-break}^2 \right)
///     \right)
/// \f]
///
/// Where:
/// - \f$B = U - I\f$ Biot strain
/// - \f$B_{sym} = (1./N_{G_1}) * \sum_i ( G_1(i) * B * G_1(i)^{-1} )\f$,
/// - \f$B_{sym-break} = B - B_{sym}\f$, the symmetry-breaking Biot strain
/// - \f$G_1(i)\f$: Element \f$i\f$ of the reference lattice point group
/// - \f$N_{G_1}\f$: Parent point group size
/// - \f$B_{reverse-sym-break} = B_{reverse} - B_{reverse-sym}\f$ is calculated
///   similarly to \f$B_{sym-break}\f$, but using \f$B_{reverse} = V_{reverse}
///   - I\f$ in place of \f$B\f$.
///
/// \param deformation_gradient The deformation gradient.
/// \param lattice1_point_group Reference lattice point group. Use the crystal
/// point group
///     of the prim if mapping structures.
///
double symmetry_breaking_strain_cost(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group) {
  Eigen::Matrix3d const &F = deformation_gradient;
  Eigen::Matrix3d U_sym_breaking =
      symmetry_breaking_right_stretch(F, lattice1_point_group);
  Eigen::Matrix3d V_reverse_sym_breaking = U_sym_breaking.inverse();

  // M.squaredNorm() = squared Frobenius norm of M = tr(M*M.transpose())
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  return ((V_reverse_sym_breaking - I).squaredNorm() +
          (U_sym_breaking - I).squaredNorm()) /
         6.;
}

}  // namespace mapping
}  // namespace CASM
