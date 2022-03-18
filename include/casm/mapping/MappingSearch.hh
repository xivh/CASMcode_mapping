#ifndef CASM_mapping_MappingSearch
#define CASM_mapping_MappingSearch

#include <map>
#include <set>

#include "casm/mapping/AtomMapping.hh"
#include "casm/mapping/LatticeMapping.hh"
#include "casm/mapping/SearchData.hh"
#include "casm/mapping/StructureMapping.hh"
#include "casm/mapping/murty.hh"

namespace CASM {
namespace mapping {

/// \brief A "node" in the search for optimal structure mappings
///
/// This encodes a particular prim, lattice, and atom mapping,
/// and includes the information needed to continue searching
/// for suboptimal assignments.
struct MappingNode {
  /// \brief Constructor
  MappingNode(
      double _lattice_cost,
      std::shared_ptr<LatticeMappingSearchData const> _lattice_mapping_data,
      double _atom_cost,
      std::shared_ptr<AtomMappingSearchData const> _atom_mapping_data,
      murty::Node _assignment_node, AtomMapping _atom_mapping,
      double _total_cost);

  /// \brief The lattice mapping cost
  double const lattice_cost;

  /// \brief Holds lattice mapping-specific data used
  ///     for mapping searches
  std::shared_ptr<LatticeMappingSearchData const> const lattice_mapping_data;

  /// \brief The atom mapping cost
  double const atom_cost;

  /// \brief Data that can be used for all atom mappings with the
  ///     same lattice mapping and trial translation
  std::shared_ptr<AtomMappingSearchData const> const atom_mapping_data;

  /// \brief Encodes a constrained solution to the atom-to-site
  ///     assignment problem
  ///
  /// This includes the information necessary (assignments forced
  /// on and forced off) to continue searching for suboptimal
  /// assignments. When solved, the sub-assignment problem is
  /// stored in assignment_node.sub_assignment.
  murty::Node const assignment_node;

  /// \brief AtomMapping solution obtained from assignment_node
  ///
  /// This includes the atom mapping permutation, displacements,
  /// and translation determined from:
  /// - the lattice mapping, stored in lattice_mapping_data
  /// - the trial translation and site displacements, stored in
  ///   atom_mapping_data
  /// - the constrained assignment problem solution stored in
  ///   assignment_node
  AtomMapping const atom_mapping;

  /// \brief The total mapping cost
  double const total_cost;

  /// \brief Compare by total_cost only
  bool operator<(MappingNode const &rhs) const {
    return this->total_cost < rhs.total_cost;
  }
};

/// \brief Performs structure mapping searches
struct MappingSearch {
  /// \brief Constructor
  MappingSearch(double _min_cost = 0.0, double _max_cost = 1e20,
                double _lattice_cost_weight = 0.5,
                std::string _atom_cost_method = "isotropic_atom_cost",
                int _k_best = 1, double _infinity = 1e20,
                double _cost_tol = 1e-5);

  /// \brief A queue of structure mappings, sorted by total
  ///     cost only
  ///
  /// This stores mappings with exactly repeated total cost
  /// in order of insertion.
  std::multiset<MappingNode> queue;

  /// \brief Results, sorted by total cost, satisifying the
  ///     min/max cost and k-best criteria
  ///
  /// This stores mappings with exactly repeated total cost
  /// in order of insertion.
  std::multimap<StructureMappingCost, StructureMapping> results;

  /// \brief Holds results that are approximately tied with
  ///     the k-best result in results
  ///
  /// This stores mappings with exactly repeated total cost
  /// in order of insertion.
  std::multimap<StructureMappingCost, StructureMapping> overflow;

  /// \brief Minimum cost result to keep
  double min_cost;

  /// \brief Maximum cost result to keep
  double max_cost;

  /// \brief Used to determine lattice_cost and components of the total cost
  ///
  /// The total cost is:
  ///
  ///     lattice_cost_weight * lattice_cost +
  ///     (1. - lattice_cost_weight) * atom_cost
  ///
  double lattice_cost_weight;

  /// \brief If true, use symmetry-break atom cost, else use isotropic atom cost
  bool is_symmetry_breaking_atom_cost;

  /// \brief Maximum number of results to keep (approximate ties are also kept)
  int k_best;

  /// \brief Cost used to prevent assignments that are not allowed
  double infinity;

  /// \brief Tolerance used for comparing costs
  double cost_tol;

  double make_atom_cost(
      Eigen::MatrixXd const &displacement,
      LatticeMappingSearchData const &lattice_mapping_data) const;

  double make_total_cost(double lattice_cost, double atom_cost) const;

  /// \brief Make a new MappingNode from an assignment problem node
  ///     with a solved sub_assignment
  MappingNode make_mapping_node_from_assignment_node(
      murty::Node assignment_node, double lattice_cost,
      std::shared_ptr<LatticeMappingSearchData const> lattice_mapping_data,
      std::shared_ptr<AtomMappingSearchData const> atom_mapping_data);

  /// \brief Insert mapping node into this->queue & this->results,
  ///     maintaining k-best results
  std::multiset<MappingNode>::iterator insert(MappingNode mapping_node);

  /// \brief Make assignment and insert mapping node
  ///     into this->queue & this->results, maintaining k-best results
  std::multiset<MappingNode>::iterator make_and_insert_mapping_node(
      double lattice_cost,
      std::shared_ptr<LatticeMappingSearchData const> lattice_mapping_data,
      std::shared_ptr<AtomMappingSearchData const> atom_mapping_data,
      std::map<Index, Index> forced_on = {},
      std::vector<std::pair<Index, Index>> forced_off = {});

  /// \brief Make the next level of sub-optimal assignments and
  ///     inserts them into this->queue & this->results, maintaining
  ///     k-best results
  std::vector<std::multiset<MappingNode>::iterator> partition(
      MappingNode const &node);
};

}  // namespace mapping
}  // namespace CASM

#endif
