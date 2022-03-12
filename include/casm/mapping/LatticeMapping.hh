#ifndef CASM_mapping_LatticeMapping
#define CASM_mapping_LatticeMapping

#include "casm/global/eigen.hh"

namespace CASM {
namespace mapping {

// Note: See source file for full documentation

/// \brief Lattice mapping transformation
struct LatticeMapping {
  /// \brief Constructor, see class documentation for definitions
  LatticeMapping(Eigen::Matrix3d const &_deformation_gradient,
                 Eigen::Matrix3d const &_transformation_matrix_to_super,
                 Eigen::Matrix3d const &_reorientation);

  /// \brief Parent-to-child deformation gradient tensor, F
  Eigen::Matrix3d deformation_gradient;

  /// \brief Parent lattice transformation_matrix_to_super, T
  Eigen::Matrix3d transformation_matrix_to_super;

  /// \brief Parent superlattice reorienation matrix, N
  Eigen::Matrix3d reorientation;

  /// \brief U, of F = Q * U = V * Q
  Eigen::Matrix3d right_stretch;

  /// \brief Rigid transformation (isometry) matrix Q, of F = Q * U = V * Q
  ///
  /// Note: Q.inverse() == Q.transpose()
  Eigen::Matrix3d isometry;

  /// \brief V, of F = Q * U = V * Q
  Eigen::Matrix3d left_stretch;
};

}  // namespace mapping
}  // namespace CASM

#endif
