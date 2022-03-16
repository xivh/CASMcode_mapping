#ifndef CASM_mapping_StructureMapping
#define CASM_mapping_StructureMapping

#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"

namespace CASM {
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

}  // namespace mapping
}  // namespace CASM

#endif
