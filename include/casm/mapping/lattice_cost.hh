#ifndef CASM_mapping_lattice_cost
#define CASM_mapping_lattice_cost

#include <vector>

#include "casm/global/eigen.hh"

namespace CASM {

namespace xtal {
struct SymOp;
}  // namespace xtal

namespace mapping {

/// \brief Returns the volume-normalized strain cost, calculated to be
/// invariant to which structure is the child/parent.
double isotropic_strain_cost(Eigen::Matrix3d const &deformation_gradient);

/// \brief Returns the symmetrized right stretch tensor
Eigen::Matrix3d symmetrized_right_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group);

/// \brief Returns the symmetry-breaking component of the right stretch tensor
Eigen::Matrix3d symmetry_breaking_right_stretch(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group);

/// \brief Returns the symmetry-breaking strain cost, calculated to be
/// invariant to which structure is the child/parent.
double symmetry_breaking_strain_cost(
    Eigen::Matrix3d const &deformation_gradient,
    std::vector<xtal::SymOp> const &lattice1_point_group);

}  // namespace mapping
}  // namespace CASM

#endif
