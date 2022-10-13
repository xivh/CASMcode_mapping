#ifndef CASM_mapping_LatticeMapping
#define CASM_mapping_LatticeMapping

#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {
namespace xtal {
class Lattice;
}

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

struct ScoredLatticeMapping : public LatticeMapping {
  ScoredLatticeMapping(double _lattice_cost, LatticeMapping _lattice_mapping)
      : LatticeMapping(_lattice_mapping), lattice_cost(_lattice_cost) {}

  double lattice_cost;
};

struct LatticeMappingResults {
  typedef std::vector<ScoredLatticeMapping>::size_type size_type;
  typedef std::vector<ScoredLatticeMapping>::const_iterator const_iterator;

  size_type size() const { return data.size(); }

  const_iterator begin() const { return data.begin(); }

  const_iterator end() const { return data.end(); }

  std::vector<ScoredLatticeMapping> data;
};

}  // namespace mapping
}  // namespace CASM

#endif
