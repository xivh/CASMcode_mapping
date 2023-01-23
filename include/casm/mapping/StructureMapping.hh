#ifndef CASM_mapping_StructureMapping
#define CASM_mapping_StructureMapping

#include <memory>
#include <vector>

#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"

namespace CASM {
namespace xtal {
class BasicStructure;
class SimpleStructure;
struct SymOp;
}  // namespace xtal
namespace mapping {

/// \brief Used to sort structure mappings and store the lattice
///     mapping cost, atom mapping cost, and combined total cost.
struct StructureMappingCost {
  StructureMappingCost(double _lattice_cost, double _atom_cost,
                       double _total_cost)
      : lattice_cost(_lattice_cost),
        atom_cost(_atom_cost),
        total_cost(_total_cost) {}
  double lattice_cost;
  double atom_cost;
  double total_cost;

  bool operator<(StructureMappingCost const &rhs) const {
    return this->total_cost < rhs.total_cost;
  }
};

/// \brief Structure mapping transformation
///
/// Combines prim, lattice mapping, and atom mapping to form a complete
/// structure mapping transformation
struct StructureMapping {
  StructureMapping(
      std::shared_ptr<xtal::BasicStructure const> const &_shared_prim,
      LatticeMapping const &_lattice_mapping, AtomMapping const &_atom_mapping)
      : shared_prim(_shared_prim),
        lattice_mapping(_lattice_mapping),
        atom_mapping(_atom_mapping) {}

  std::shared_ptr<xtal::BasicStructure const> shared_prim;
  LatticeMapping lattice_mapping;
  AtomMapping atom_mapping;
};

/// \brief Return mappings that result in structures along the
///     transformation pathway from the parent to the aligned child
///     structure
StructureMapping interpolated_mapping(StructureMapping const &structure_mapping,
                                      double interpolation_factor);

struct ScoredStructureMapping : public StructureMapping {
  ScoredStructureMapping(double _lattice_cost, double _atom_cost,
                         double _total_cost,
                         StructureMapping _structure_mapping)
      : StructureMapping(_structure_mapping),
        lattice_cost(_lattice_cost),
        atom_cost(_atom_cost),
        total_cost(_total_cost) {}

  ScoredStructureMapping(StructureMappingCost const &_cost,
                         StructureMapping _structure_mapping)
      : ScoredStructureMapping(_cost.lattice_cost, _cost.atom_cost,
                               _cost.total_cost, _structure_mapping) {}

  ScoredStructureMapping(
      std::pair<StructureMappingCost, StructureMapping> const &_pair)
      : ScoredStructureMapping(_pair.first, _pair.second) {}

  double lattice_cost;
  double atom_cost;
  double total_cost;
};

struct StructureMappingResults {
  typedef std::vector<ScoredStructureMapping>::size_type size_type;
  typedef std::vector<ScoredStructureMapping>::const_iterator const_iterator;

  StructureMappingResults(std::vector<ScoredStructureMapping> _data = {})
      : data(_data) {}

  size_type size() const { return data.size(); }

  const_iterator begin() const { return data.begin(); }

  const_iterator end() const { return data.end(); }

  std::vector<ScoredStructureMapping> data;
};

/// \brief Return the mapped structure, with implied vacancies,
///     strain, and atomic displacement
xtal::SimpleStructure make_mapped_structure(
    StructureMapping const &structure_mapping,
    xtal::SimpleStructure const &unmapped_structure);

/// \brief Make the structure mapping to a different, but equivalent
///     structure
StructureMapping make_mapping_to_equivalent_structure(
    xtal::SymOp const &op, xtal::Lattice const &target,
    StructureMapping const &structure_mapping);

/// \brief Make an equivalent mapping to a different, but symmetrically
///     equivalent superlattice of the prim
StructureMapping make_mapping_to_equivalent_superlattice(
    xtal::Lattice const &target, StructureMapping const &structure_mapping,
    std::vector<xtal::SymOp> const &group);

}  // namespace mapping
}  // namespace CASM

#endif
