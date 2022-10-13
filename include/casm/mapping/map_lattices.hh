#ifndef CASM_mapping_map_lattices
#define CASM_mapping_map_lattices

#include <optional>
#include <string>
#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class Lattice;
struct SymOp;
}  // namespace xtal

namespace mapping {
struct LatticeMapping;
struct LatticeMappingResults;

// Note: See source file for full documentation

/// \brief Find lattice mappings
LatticeMappingResults map_lattices(
    xtal::Lattice const &lattice1, xtal::Lattice const &lattice2,
    std::optional<Eigen::Matrix3d> T = std::nullopt,
    int reorientation_range = 1,
    std::vector<xtal::SymOp> lattice1_point_group = std::vector<xtal::SymOp>{},
    std::vector<xtal::SymOp> lattice2_point_group = std::vector<xtal::SymOp>{},
    double min_cost = 0.0, double max_cost = 1e20,
    std::string cost_method = std::string("isotropic_strain_cost"),
    std::optional<int> k_best = std::nullopt, double cost_tol = 1e-5);

}  // namespace mapping
}  // namespace CASM

#endif
