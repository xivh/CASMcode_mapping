#ifndef CASM_mapping_atom_cost
#define CASM_mapping_atom_cost

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
class BasicStructure;
class UnitCellCoordIndexConverter;
}  // namespace xtal

namespace mapping {

struct LatticeMapping;

/// \brief A volume-normalized mean-squared-displacement cost
double make_geometric_atom_cost(
    Eigen::Matrix3d const &supercell_lattice_column_vector_matrix,
    Eigen::MatrixXd const &displacement);

/// \brief Return the "isotropic" atom cost
double make_isotropic_atom_cost(
    Eigen::Matrix3d const &prim_lattice_column_vector_matrix,
    LatticeMapping const &lattice_mapping, Eigen::MatrixXd const &displacement);

/// \brief Return the symmetry-preserving component of displacement
Eigen::MatrixXd make_symmetry_preserving_displacement(
    Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes);

/// \brief Return the symmetry-breaking component of displacement
Eigen::MatrixXd make_symmetry_breaking_displacement(
    Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes);

/// \brief Return the "symmetry breaking atom cost"
double make_symmetry_breaking_atom_cost(
    Eigen::Matrix3d const &prim_lattice_column_vector_matrix,
    LatticeMapping const &lattice_mapping, Eigen::MatrixXd const &displacement,
    xtal::UnitCellCoordIndexConverter const &unitcellcoord_index_converter,
    std::vector<Eigen::MatrixXd> const &prim_sym_invariant_displacement_modes);

}  // namespace mapping
}  // namespace CASM

#endif
