#ifndef CASM_mapping_StructureMapping
#define CASM_mapping_StructureMapping

#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"

namespace CASM {
namespace mapping {

/// \brief Structure mapping transformation
///
/// Combines lattice mapping and atom mapping to form a complete structure
/// mapping transformation
struct StructureMapping {
  StructureMapping(LatticeMapping const &_lattice_mapping,
                   AtomMapping const &_atom_mapping)
      : lattice_mapping(_lattice_mapping), atom_mapping(_atom_mapping) {}

  LatticeMapping lattice_mapping;
  AtomMapping atom_mapping;
};

}  // namespace mapping
}  // namespace CASM

#endif
