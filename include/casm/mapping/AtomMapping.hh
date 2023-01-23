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

/// \brief Return mappings that result in atom positions along the
///     transformation pathway from the parent to the aligned child
///     structure
AtomMapping interpolated_mapping(AtomMapping const &atom_mapping,
                                 double interpolation_factor);

struct ScoredAtomMapping : public AtomMapping {
  ScoredAtomMapping(double _atom_cost, AtomMapping _atom_mapping)
      : AtomMapping(_atom_mapping), atom_cost(_atom_cost) {}

  double atom_cost;
};

struct AtomMappingResults {
  typedef std::vector<ScoredAtomMapping>::size_type size_type;
  typedef std::vector<ScoredAtomMapping>::const_iterator const_iterator;

  AtomMappingResults(std::vector<ScoredAtomMapping> _data = {}) : data(_data) {}

  size_type size() const { return data.size(); }

  const_iterator begin() const { return data.begin(); }

  const_iterator end() const { return data.end(); }

  std::vector<ScoredAtomMapping> data;
};

}  // namespace mapping
}  // namespace CASM

#endif
