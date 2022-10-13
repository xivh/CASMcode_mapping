#ifndef CASM_mapping_map_atoms
#define CASM_mapping_map_atoms

#include <optional>
#include <string>
#include <vector>

namespace CASM {

namespace xtal {
class BasicStructure;
class Lattice;
class SimpleStructure;
struct SymOp;
}  // namespace xtal

namespace mapping {
struct AtomMapping;
struct AtomMappingResults;
struct LatticeMapping;

// Note: See source file for full documentation

/// \brief Find atom mappings
///
/// \param prim The reference "parent" structure
/// \param structure2 The "child" structure
/// \param lattice_mapping The lattice mapping to use
/// \param prim_factor_group Used to skip mappings that are
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param min_cost Keep results with total cost >= min_cost
/// \param max_cost Keep results with total cost <= max_cost
/// \param atom_cost_method One of "isotropic_atom_cost" or
///     "symmetry_breaking_atom_cost"
/// \param k_best Keep the k_best results satisfying the min_cost and
///     max_cost constraints. If there are approximate ties, those
///     will also be kept.
/// \param cost_tol Tolerance for checking if lattice mapping costs are
///     approximately equal
AtomMappingResults map_atoms(
    xtal::BasicStructure const &prim, xtal::SimpleStructure const &structure2,
    LatticeMapping const &lattice_mapping,
    std::vector<xtal::SymOp> prim_factor_group = std::vector<xtal::SymOp>{},
    double min_cost = 0.0, double max_cost = 1e20,
    std::string atom_cost_method = std::string("isotropic_atom_cost"),
    int k_best = 1, double cost_tol = 1e-5);

}  // namespace mapping
}  // namespace CASM

#endif
