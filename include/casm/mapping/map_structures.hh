#ifndef CASM_mapping_map_structures
#define CASM_mapping_map_structures

#include <optional>
#include <string>
#include <vector>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class SimpleStructure;
struct SymOp;
}  // namespace xtal

namespace mapping {
struct StructureMapping;
struct StructureMappingResults;

// Note: See source file for full documentation

/// \brief Find structure mappings, given a range of parent superstructure
/// volumes
StructureMappingResults map_structures(
    xtal::BasicStructure const &prim, xtal::SimpleStructure const &structure2,
    Index max_vol,
    std::vector<xtal::SymOp> prim_factor_group = std::vector<xtal::SymOp>{},
    std::vector<xtal::SymOp> structure2_factor_group =
        std::vector<xtal::SymOp>{},
    Index min_vol = 1, double min_cost = 0.0, double max_cost = 1e20,
    double lattice_cost_weight = 0.5,
    std::string lattice_cost_method = std::string("isotropic_strain_cost"),
    std::string atom_cost_method = std::string("isotropic_disp_cost"),
    int k_best = 1, double cost_tol = 1e-5);

/// \brief Find structure mappings, given a range of parent superstructure
/// volumes
StructureMappingResults map_structures_v2(
    xtal::BasicStructure const &prim, xtal::SimpleStructure const &structure2,
    Index max_vol,
    std::vector<xtal::SymOp> prim_factor_group = std::vector<xtal::SymOp>{},
    std::vector<xtal::SymOp> structure2_factor_group =
        std::vector<xtal::SymOp>{},
    Index min_vol = 1, double min_cost = 0.0, double max_cost = 1e20,
    double lattice_cost_weight = 0.5,
    std::string lattice_cost_method = std::string("isotropic_strain_cost"),
    std::string atom_cost_method = std::string("isotropic_disp_cost"),
    int k_best = 1, double cost_tol = 1e-5);

}  // namespace mapping
}  // namespace CASM

#endif
