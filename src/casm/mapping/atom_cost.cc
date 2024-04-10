#include "casm/mapping/atom_cost.hh"

#include "casm/crystallography/LinearIndexConverter.hh"
#include "casm/global/eigen.hh"
#include "casm/mapping/LatticeMapping.hh"

namespace CASM {
namespace mapping {

/// \brief A volume-normalized mean-squared-displacement cost
///
/// The geometric atom cost is:
/// \f[
///     \left( \frac{3v}{4\pi} \right)^{-2/3} \frac{\sum^N_i d_{i}^2}{N}
/// \f]
/// where:
///
/// - d_{i}: The site-to-atom displacement, F*(r1[i] + d[i]) = r2[perm[i]] +
/// trans
/// - v: The volume per site
///
/// \param supercell_lattice_column_vector_matrix A shape=(3,3) matrix
///     containing the supercell lattice vectors as columns.
/// \param displacement Shape=(3, N_supercell_site) matrix with
///     the site-to-atom displacements, as defined in an
///     AtomMapping.
double make_geometric_atom_cost(
    Eigen::Matrix3d const &supercell_lattice_column_vector_matrix,
    Eigen::MatrixXd const &displacement) {
  double N_site = displacement.cols();
  double volume_per_site =
      std::abs(supercell_lattice_column_vector_matrix.determinant()) / N_site;
  return std::pow(3. * volume_per_site / (4. * M_PI), -2. / 3.) *
         displacement.squaredNorm() / N_site;
}

/// \brief Return the "isotropic" atom cost
///
///
/// Want to get an atom cost that:
/// - does not depend on which structure is "parent" and
///   which is "child"
/// - is volume-normalized
///
/// Determine "parent" lattice S1, and "child" lattice, L2, from
/// lattice_mapping:
///
///     // U * L1 * T * N == L2
///     // U * S1 == L2
///     // -->
///     S1 = L1 * T * N
///     L2 = U * S1
///
/// Determine "parent" displacements, d, and "child" displacements,
/// d_reverse, from the lattice_mapping and atom_mapping:
///
///     // U (r1[i] + d) = r2[i] + trans
///     // r1 + trans_reverse = U_reverse * (r2 + d_reverse)
///     // -->
///     U_reverse = U.inv
///     d_reverse = -U * d
///     trans_reverse = -U.inv * trans
///
/// This returns:
///
///     make_geometric_atom_cost(S1, d)*0.5 + make_geometric_atom_cost(L2,
///     d_reverse)
///
/// \param prim_lattice_column_vector_matrix A shape=(3,3) matrix
///     containing the prim lattice vectors as columns.
/// \param lattice_mapping A LatticeMappng solution
/// \param displacement Shape=(3, N_supercell_site) matrix with
///     the site-to-atom displacements, as defined in an
///     AtomMapping.
///
/// \return The isotropic atom cost
///
double make_isotropic_atom_cost(
    Eigen::Matrix3d const &prim_lattice_column_vector_matrix,
    LatticeMapping const &lattice_mapping,
    Eigen::MatrixXd const &displacement) {
  Eigen::Matrix3d const &L1 = prim_lattice_column_vector_matrix;
  Eigen::Matrix3d const &T = lattice_mapping.transformation_matrix_to_super;
  Eigen::Matrix3d const &N = lattice_mapping.reorientation;
  Eigen::Matrix3d const &U = lattice_mapping.right_stretch;

  Eigen::Matrix3d S1 = L1 * T * N;
  Eigen::Matrix3d L2 = U * S1;
  Eigen::MatrixXd const &d = displacement;
  Eigen::MatrixXd d_reverse = -U * d;

  double isotropic_atom_cost = (make_geometric_atom_cost(S1, d) +
                                make_geometric_atom_cost(L2, d_reverse)) /
                               2.;
  return isotropic_atom_cost;
}

/// \brief Return the symmetry-preserving component of displacement
///
/// \param displacement Shape=(3, N_supercell_site) matrix with
///     the site-to-atom displacements, as defined in an
///     AtomMapping.
/// \param unitcellcoord_index_converter Gives the coordinates
///     for sites associated with the columns of the supercell
///     displacement matrix.
/// \param prim_sym_invariant_displacement_modes Size=N_mode, with
///     shape=(3,N_prim_site) matrices, giving the
///     symmetry invariant displacement modes, with columns
///     giving the displacement associated with each prim site
///     for a given mode.
Eigen::MatrixXd make_symmetry_preserving_displacement(
    Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes) {
  Eigen::Map<Eigen::VectorXd const> unrolled_displacement(displacement.data(),
                                                          displacement.size());
  Index N_site = displacement.cols();

  Eigen::MatrixXd symmetry_preserving_displacement =
      Eigen::MatrixXd::Zero(3, displacement.cols());

  for (auto const &prim_mode : prim_sym_invariant_displacement_modes) {
    // Make supercell mode corresponding to prim mode
    Eigen::MatrixXd supercell_mode = Eigen::MatrixXd::Zero(3, N_site);
    for (Index l = 0; l < N_site; ++l) {
      Index b = unitcellcoord_index_converter(l).sublattice();
      supercell_mode.col(l) = prim_mode.col(b);
    }
    supercell_mode /= supercell_mode.norm();

    // Project the displacement mode on to this symmetry invariant mode
    Eigen::VectorXd unrolled_supercell_mode = Eigen::Map<Eigen::VectorXd const>(
        supercell_mode.data(), supercell_mode.size());
    double projected_comp = unrolled_displacement.dot(unrolled_supercell_mode);
    symmetry_preserving_displacement += projected_comp * supercell_mode;
  }
  return symmetry_preserving_displacement;
}

/// \brief Return the symmetry-breaking component of displacement
///
/// \param displacement Shape=(3, N_supercell_site) matrix with
///     the site-to-atom displacements, as defined in an
///     AtomMapping.
/// \param unitcellcoord_index_converter Gives the coordinates
///     for sites associated with the columns of the supercell
///     displacement matrix.
/// \param prim_sym_invariant_displacement_modes Size=N_mode, with
///     shape=(3,N_prim_site) matrices, giving the
///     symmetry invariant displacement modes, with columns
///     giving the displacement associated with each prim site
///     for a given mode.
Eigen::MatrixXd make_symmetry_breaking_displacement(
    Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes) {
  return displacement - make_symmetry_preserving_displacement(
                            displacement, unitcellcoord_index_converter,
                            prim_sym_invariant_displacement_modes);
}

/// \brief Return the "symmetry breaking atom cost"
///
/// This is equal to the isotropic atom cost due to the
/// symmetry-breaking displacements.
///
/// \param prim_lattice_column_vector_matrix A shape=(3,3) matrix
///     containing the prim lattice vectors as columns.
/// \param lattice_mapping A LatticeMappng solution
/// \param displacement Shape=(3, N_supercell_site) matrix with
///     the site-to-atom displacements, as defined in an
///     AtomMapping.
/// \param unitcellcoord_index_converter Gives the coordinates
///     for sites associated with the columns of the supercell
///     displacement matrix.
/// \param prim_sym_invariant_displacement_modes Size=N_mode, with
///     shape=(3,N_prim_site) matrices, giving the
///     symmetry invariant displacement modes, with columns
///     giving the displacement associated with each prim site
///     for a given mode.
double make_symmetry_breaking_atom_cost(
    Eigen::Matrix3d const &prim_lattice_column_vector_matrix,
    LatticeMapping const &lattice_mapping, Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes) {
  Eigen::MatrixXd symmetry_breaking_displacement =
      make_symmetry_breaking_displacement(
          displacement, unitcellcoord_index_converter,
          prim_sym_invariant_displacement_modes);

  return make_isotropic_atom_cost(prim_lattice_column_vector_matrix,
                                  lattice_mapping,
                                  symmetry_breaking_displacement);
}

}  // namespace mapping
}  // namespace CASM
