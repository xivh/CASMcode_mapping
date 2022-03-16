#include "casm/mapping/LatticeMapping.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/misc/CASM_Eigen_math.hh"

namespace CASM {
namespace mapping {

/// \struct LatticeMapping
/// \brief Lattice mapping transformation
///
/// The lattice mapping transformation has the form:
///
/// \f[
///     F L_1 T N = L_2.
/// \f]
///
/// where \f$L_1\f$ are the reference "parent" lattice vectors, and
/// \f$L_2\f$ are the "child" lattice vectors being mapped to the
/// parent vectors, as columns of shape=(3,3) matrices. The other
/// shape=(3,3) matrices are:
///
/// - \f$F\f$, the shape=(3,3) parent-to-child deformation
///   gradient tensor
/// - \f$T\f$, an integer transformation matrix that generates a
///   superlattice of \f$L_1\f$
/// - \f$N\f$, a unimodular reorientation matrix that generates a
///   lattice equivalent to \f$L_1 T\f$ with reoriented lattice
///   vectors
///
/// There are infinitely many potential choices of N, which leads to
/// infinitely many potential lattice mapping transformations between
/// two lattices. A cost function of the stretch, \f$V\f$, is used
/// to score mappings.
///

/// \brief Constructor
///
/// \param _deformation_gradient Parent-to-child deformation gradient
///     tensor, F
/// \param _transformation_matrix_to_super Parent lattice
///     transformation_matrix_to_super, T
/// \param _reorientation Parent superlattice reorienation matrix, N.
///     Must be unimodular (approximately integer valued, with
///     determinant = +/-1).
///
/// \throws If _reorientation is not unimodular
///
LatticeMapping::LatticeMapping(
    Eigen::Matrix3d const &_deformation_gradient,
    Eigen::Matrix3d const &_transformation_matrix_to_super,
    Eigen::Matrix3d const &_reorientation)
    : deformation_gradient(_deformation_gradient),
      transformation_matrix_to_super(_transformation_matrix_to_super),
      reorientation(_reorientation),
      right_stretch(polar_decomposition(deformation_gradient)),
      isometry(deformation_gradient * right_stretch.inverse()),
      left_stretch(deformation_gradient * isometry.transpose()) {
  if (!is_unimodular(reorientation, 1e-5)) {
    throw std::runtime_error(
        "Error in LatticeMapping: reorientation matrix is not unimodular");
  }
}

/// \brief Return mapped lattice, L2 = F * L1 * T * N
xtal::Lattice make_mapped_lattice(xtal::Lattice const &prim_lattice,
                                  LatticeMapping const &lattice_mapping) {
  auto const &F = lattice_mapping.deformation_gradient;
  auto const &L1 = prim_lattice.lat_column_mat();
  auto const &T = lattice_mapping.transformation_matrix_to_super;
  auto const &N = lattice_mapping.reorientation;
  return xtal::Lattice(F * L1 * T * N, prim_lattice.tol());
}

}  // namespace mapping
}  // namespace CASM
