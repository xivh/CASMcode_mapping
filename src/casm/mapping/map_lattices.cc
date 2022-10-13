#include "casm/mapping/map_lattices.hh"

#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/SymType.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/impl/LatticeMap.hh"
#include "casm/mapping/misc.hh"

namespace CASM {
namespace mapping_impl {

/// \brief Used to sort lattice mappings by cost, then CASM lattice comparison
struct LatticeMappingKey {
  LatticeMappingKey(double _cost, double _cost_tol, xtal::Lattice _lattice)
      : cost(_cost), cost_tol(_cost_tol), lattice(_lattice) {}

  double cost;
  double cost_tol;
  xtal::Lattice lattice;

  bool operator<(LatticeMappingKey const &RHS) const {
    if (this->cost < RHS.cost - cost_tol) {
      return true;
    }
    if (this->cost > RHS.cost + cost_tol) {
      return false;
    }
    return this->lattice > RHS.lattice;
  }
};

}  // namespace mapping_impl

namespace mapping {

/// \brief Find lattice mappings
///
/// Find and rank lattice mappings (see `LatticeMapping`) by searching
/// over equivalent lattices with different lattice vectors by varying
/// the reorientation matrix, N.
///
/// This method is often used inside a loop over superlattice
/// transformation matrices, T, to check mappings between distinct
/// superlattices.
///
/// For strain cost definitions, see Python documentation.
///
/// For more details, see J.C. Thomas, A.R. Natarajan, and A.V. Van der Ven,
/// npj Computational Materials (2021)7:164;
/// https://doi.org/10.1038/s41524-021-00627-0
///
/// \param lattice1 The referece "parent" lattice
/// \param lattice2 The "child" lattice
/// \param T Mapping is performed for L1 * T * N <-> L2, where T is an
///     approximately integer transformation matrix to a supercell of L1.
///     The default value is the identity matrix.
/// \param reorientation_range The absolute value of the maximum element in
///     reorientation matrix, N. This determines how many equivalent lattice
///     vector reorientations are checked. Usually 1 is sufficient.
/// \param L1_point_group Used to skip reorientation matrices that result in
///     symmetrically equivalent mappings. The default (empty), is
///     equivalent to only including the identity operation.
/// \param L2_point_group Used to skip reorientation matrices that
///     result in symmetrically equivalent mappings. The default (empty), is
///     equivalent to just including the identity operation.
/// \param min_cost Keep results with cost >= min_cost
/// \param max_cost Keep results with cost <= max_cost
/// \param cost_method One of "isotropic_strain_cost" or
///     "symmetry_breaking_strain_cost"
/// \param k_best If k_best.has_value(), then only keep the k_best results
///     satisfying the min_cost and max_cost constraints. If there are
///     approximate ties, those will also be kept.
/// \param cost_tol Tolerance for checking if lattice mapping costs are
///     approximately equal
///
/// \returns A vector of {cost, lattice_mapping}, giving lattice
///     mapping solutions and their costs, sorted by lattice mapping cost.
///
LatticeMappingResults map_lattices(
    xtal::Lattice const &lattice1, xtal::Lattice const &lattice2,
    std::optional<Eigen::Matrix3d> T, int reorientation_range,
    std::vector<xtal::SymOp> lattice1_point_group,
    std::vector<xtal::SymOp> lattice2_point_group, double min_cost,
    double max_cost, std::string cost_method, std::optional<int> k_best,
    double cost_tol) {
  double init_better_than = 1e20;
  bool symmetrize_strain_cost;
  if (cost_method == "isotropic_strain_cost") {
    symmetrize_strain_cost = false;
  } else if (cost_method == "symmetry_breaking_strain_cost") {
    symmetrize_strain_cost = true;
  } else {
    throw std::runtime_error(
        "Error in map_lattices: cost_method not recognized");
  }
  if (k_best.has_value() && *k_best < 1) {
    throw std::runtime_error(
        "Error in map_lattices: k_best < 1 is not allowed");
  }
  if (lattice1_point_group.empty()) {
    lattice1_point_group.push_back(xtal::SymOp::identity());
  }
  if (lattice2_point_group.empty()) {
    lattice2_point_group.push_back(xtal::SymOp::identity());
  }
  if (!T.has_value()) {
    T = Eigen::Matrix3d::Identity();
  }
  Eigen::Matrix3d L1 = lattice1.lat_column_mat();
  xtal::Lattice parent_superlattice(L1 * T.value(), lattice1.tol());
  mapping_impl::LatticeMap latmap(
      parent_superlattice, lattice2, reorientation_range, lattice1_point_group,
      lattice2_point_group, init_better_than, symmetrize_strain_cost, cost_tol);

  // the k-best results
  std::multimap<mapping_impl::LatticeMappingKey, LatticeMapping> results;

  // results that are approximately equal to last-place result
  std::multimap<mapping_impl::LatticeMappingKey, LatticeMapping> overflow;

  while (latmap) {
    double cost = latmap.strain_cost();
    if (cost > (min_cost - cost_tol) && cost < (max_cost + cost_tol)) {
      results.emplace(
          mapping_impl::LatticeMappingKey(
              cost, cost_tol, xtal::Lattice(L1 * T.value() * latmap.matrixN())),
          LatticeMapping(latmap.deformation_gradient(), T.value(),
                         latmap.matrixN()));

      // maintain results.size() <= *k_best, keep approximately equal results,
      // shrinks max_cost to results.rbegin() if results.size() == *k_best
      maintain_k_best_results(
          k_best, cost_tol, results, overflow,
          [](mapping_impl::LatticeMappingKey const &key) { return key.cost; });
      if (k_best.has_value() && results.size() == *k_best) {
        max_cost = results.rbegin()->first.cost;
      }
    }

    // iterate reorientation matrices until finding one with cost < max_cost +
    // cost_tol
    latmap.next_mapping_better_than(max_cost);
  }

  while (overflow.size()) {
    results.insert(overflow.extract(overflow.begin()));
  }
  LatticeMappingResults final;
  for (auto const &pair : results) {
    final.data.emplace_back(pair.first.cost, pair.second);
  }
  return final;
}

}  // namespace mapping
}  // namespace CASM
