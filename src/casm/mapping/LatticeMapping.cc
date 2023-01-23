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

/// \brief Return mappings that result in lattices along the
///     transformation pathway from the parent to the aligned child
///     lattice
///
/// \param structure_mapping Original lattice mapping
/// \param interpolation_factor Interpolation factor. The value 0.0
///     corresponds to the child lattice mapped to the ideal parent
///     lattice; and the value 1.0 corresponds to the child
///     lattice, transformed by isometry to align with the ideal
///     parent lattice.
LatticeMapping interpolated_mapping(LatticeMapping const &lattice_mapping,
                                    double interpolation_factor) {
  double f = interpolation_factor;
  LatticeMapping const &lmap = lattice_mapping;

  auto const &Q = lmap.isometry;
  auto const &U = lmap.right_stretch;
  Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
  Eigen::Matrix3d interpolated_U = I + (U - I) * f;
  return LatticeMapping(Q * interpolated_U, lmap.transformation_matrix_to_super,
                        lmap.reorientation);
}

/// \brief Return the mapped lattice
///
/// The "mapped lattice" is constructed from the parent lattice
/// by applying the lattice mapping without isometry. It has
/// the lattice vector column matrix:
///
/// \f[
///     L_m = U L_1 T N
/// \f]
///
/// where \f$L_1\f$ is the reference "parent" lattice vectors, and
/// \f$L_m\f$ are the "mapped" lattice vectors, as columns of
/// shape=(3,3) matrices. The othe shape=(3,3) matrices are:
///
/// - \f$U\f$, the shape=(3,3) right stretch tensor of the
///   parent-to-child deformation gradient tensor
/// - \f$T\f$, an integer transformation matrix that generates a
///   superlattice of \f$L_1\f$
/// - \f$N\f$, a unimodular reorientation matrix that generates a
///   lattice equivalent to \f$L_1 T\f$ with reoriented lattice
///   vectors
///
/// \param parent_lattice The reference lattice
/// \param lattice_mapping The lattice mapping transformation
///
/// \returns mapped_lattice The mapped lattice that results
///     from transforming the parent lattice according to the
///     lattice mapping
///
xtal::Lattice make_mapped_lattice(xtal::Lattice const &parent_lattice,
                                  LatticeMapping const &lattice_mapping) {
  auto const &U = lattice_mapping.right_stretch;
  auto const &L1 = parent_lattice.lat_column_mat();
  auto const &T = lattice_mapping.transformation_matrix_to_super;
  auto const &N = lattice_mapping.reorientation;

  return xtal::Lattice(U * L1 * T * N);
}

}  // namespace mapping
}  // namespace CASM
