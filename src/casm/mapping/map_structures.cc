#include "casm/mapping/map_structures.hh"

#include "casm/crystallography/BasicStructure.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SimpleStructure.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/StructureMapping.hh"
#include "casm/mapping/impl/SimpleStrucMapCalculator.hh"
#include "casm/mapping/impl/StrucMapping.hh"

namespace CASM {
namespace mapping_impl {

// Declarations of functions defined in map_atoms.cc:

mapping::LatticeMapping make_lattice_mapping(LatticeNode const &lattice_node);

mapping::AtomMapping make_atom_mapping(
    Eigen::Matrix3d const &deformation_gradient,
    MappingNode const &mapping_node);

}  // namespace mapping_impl

namespace mapping {

/// \brief Find structure mappings, given a range of parent superstructure
/// volumes
///
/// This method finds mappings from a superstructure of a reference "parent"
/// structure to a "child" structure. It works by finding lattice mappings
/// (see `LatticeMapping`) for symmetrically unique superlattices of the
/// "parent" lattice for a range of supercell volumes, and for each
/// potential lattice mapping finding atom mappings (see `AtomMapping`).
///
/// The total structure mapping cost, total_cost, is a weighted mixture of
/// the lattice mapping cost, lattice_cost, and the atom mapping cost,
/// atom_cost:
///
///     total_cost = lattice_cost_weight*lattice_cost
///                  + (1.0 - lattice_cost_weight)*atom_cost
///
/// where lattice_cost_weight is an input parameter.
///
/// For strain and atom cost definitions, see Python documentation.
///
/// For more details, see J.C. Thomas, A.R. Natarajan, and A.V. Van der Ven,
/// npj Computational Materials (2021)7:164;
/// https://doi.org/10.1038/s41524-021-00627-0
///
/// \param prim The reference "parent" structure
/// \param structure2 The "child" structure
/// \param max_vol The maximum parent superstructure volume to consider, as
///     a multiple of the parent structure volume
/// \param prim_factor_group Used to skip mappings that are
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param structure2_factor_group Used to skip mappings that are
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param min_vol The minimum parent superstructure volume to consider, as
///     a multiple of the parent structure volume. Default=1.
/// \param min_cost Keep results with total cost >= min_cost
/// \param max_cost Keep results with total cost <= max_cost
/// \param lattice_cost_weight The fraction of the total cost due to
///     the lattice strain cost. The remaining fraction
///     (1.-lattice_cost_weight) is due to the atom cost. Default=0.5.
/// \param lattice_cost_method One of "isotropic_strain_cost" or
///     "symmetry_breaking_strain_cost"
/// \param atom_cost_method One of "isotropic_disp_cost" or
///     "symmetry_breaking_disp_cost"
/// \param k_best Keep the k_best results satisfying the min_cost and
///     max_cost constraints. If there are approximate ties, those
///     will also be kept.
/// \param cost_tol Tolerance for checking if lattice mapping costs are
///     approximately equal
StructureMappingResults map_structures(
    xtal::BasicStructure const &prim, xtal::SimpleStructure const &structure2,
    Index max_vol, std::vector<xtal::SymOp> prim_factor_group,
    std::vector<xtal::SymOp> structure2_factor_group, Index min_vol,
    double min_cost, double max_cost, double lattice_cost_weight,
    std::string lattice_cost_method, std::string atom_cost_method, int k_best,
    double cost_tol) {
  auto shared_prim = std::make_shared<xtal::BasicStructure const>(prim);

  bool symmetrize_lattice_cost;
  if (lattice_cost_method == "isotropic_strain_cost") {
    symmetrize_lattice_cost = false;
  } else if (lattice_cost_method == "symmetry_breaking_strain_cost") {
    symmetrize_lattice_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_structures: lattice_cost_method not recognized");
  }

  bool symmetrize_atom_cost;
  if (atom_cost_method == "isotropic_disp_cost") {
    symmetrize_atom_cost = false;
  } else if (atom_cost_method == "symmetry_breaking_disp_cost") {
    symmetrize_atom_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_structures: atom_cost_method not recognized");
  }

  if (k_best < 1) {
    throw std::runtime_error(
        "Error in map_structures: k_best < 1 is not allowed");
  }
  if (prim_factor_group.empty()) {
    prim_factor_group.push_back(xtal::SymOp::identity());
  }
  if (structure2_factor_group.empty()) {
    structure2_factor_group.push_back(xtal::SymOp::identity());
  }
  if (min_vol < 1) {
    throw std::runtime_error("Error in map_structures: min_vol < 1");
  }
  if (max_vol < min_vol) {
    throw std::runtime_error("Error in map_structures: max_vol < min_vol");
  }

  /// For the StrucMapper::map_deformed_struc_impose_lattice_vols method:
  /// - If invalid values of `min_vol` or `max_vol` are provided (negative
  /// values
  ///   or max_vol < min_vol), then this method throws.
  /// - Parameters `min_va_frac`, `max_va_frac`, `max_volume_change`,
  ///   and `soft_va_limit` have no effect.
  /// - `robust` search is used anyway if k_best > 1

  mapping_impl::SimpleStrucMapCalculator calculator(
      xtal::make_simple_structure(prim), prim_factor_group,
      CASM::xtal::SimpleStructure::SpeciesMode::ATOM,
      xtal::allowed_molecule_names(prim));
  double _max_volume_change = 0.5;  // no effect
  bool _robust = true;              // no effect if k_best > 1
  bool _soft_va_limit = false;      // no effect
  double _min_va_frac = 0.;         // no effect
  double _max_va_frac = 1.;         // no effect
  mapping_impl::StrucMapper strucmap(
      calculator, lattice_cost_weight, _max_volume_change, _robust,
      _soft_va_limit, cost_tol, _min_va_frac, _max_va_frac);

  if (symmetrize_lattice_cost) {
    strucmap.set_symmetrize_lattice_cost(true);
  }
  if (symmetrize_atom_cost) {
    auto prim_permute_group =
        xtal::make_permutation_representation(prim, prim_factor_group);
    strucmap.set_symmetrize_atomic_cost(true, prim_factor_group,
                                        prim_permute_group);
  }

  bool keep_invalid = false;
  std::set<mapping_impl::MappingNode> mappings =
      strucmap.map_deformed_struc_impose_lattice_vols(
          structure2, min_vol, max_vol, k_best, max_cost, min_cost,
          keep_invalid, structure2_factor_group);

  // Convert mapping_impl::MappingNode results to StructureMapping results
  StructureMappingResults results;
  for (auto const &mapping_node : mappings) {
    if (!(mapping_node.cost > (min_cost - cost_tol) &&
          mapping_node.cost < (max_cost + cost_tol))) {
      continue;
    }

    // Get LatticeMapping data
    LatticeMapping lattice_mapping =
        mapping_impl::make_lattice_mapping(mapping_node.lattice_node);

    // Get AtomMapping data
    AtomMapping atom_mapping = mapping_impl::make_atom_mapping(
        lattice_mapping.deformation_gradient, mapping_node);

    results.data.emplace_back(
        mapping_node.lattice_node.cost, mapping_node.atomic_node.cost,
        mapping_node.cost,
        StructureMapping(shared_prim, lattice_mapping, atom_mapping));
  }

  return results;
}

}  // namespace mapping
}  // namespace CASM
