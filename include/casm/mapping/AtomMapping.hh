#ifndef CASM_mapping_AtomMapping
#define CASM_mapping_AtomMapping

#include <vector>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace mapping {

// Note: See source file for full documentation

/// \brief Atom mapping transformation
struct AtomMapping {
  AtomMapping(Eigen::MatrixXd const &_displacement,
              std::vector<Index> const &_permutation,
              Eigen::Vector3d const &_translation);

  Eigen::MatrixXd displacement;
  std::vector<Index> permutation;
  Eigen::Vector3d translation;
};

}  // namespace mapping
}  // namespace CASM

#endif
