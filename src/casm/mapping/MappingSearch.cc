#include "casm/mapping/MappingSearch.hh"

#include "casm/mapping/atom_cost.hh"
#include "casm/mapping/hungarian.hh"
#include "casm/mapping/misc.hh"

namespace CASM {
namespace mapping {

namespace mapping_impl {

bool make_is_symmetry_breaking_atom_cost(std::string atom_cost_method) {
  if (atom_cost_method == "isotropic_atom_cost") {
    return false;
  } else if (atom_cost_method == "symmetry_breaking_atom_cost") {
    return true;
  } else {
    throw std::runtime_error("Error: atom_cost_method not recognized");
  }
}

/// \brief Construct an AtomMapping from an assignment solution
///
/// The site displacements are the minimum length displacements
/// that satisfy:
///
///     site_coordinate_cart[i] + site_displacements[i][j] =
///         F^{-1}*atom_coordinate_cart[j] + trial_translation
///
/// under periodic boundary conditions, and in the context of
/// a particular lattice mapping. The AtomMapping transformation
/// is:
///
///     F * (r1[i] + d[i]) = r2[perm[i]] + translation
///
/// in the context of a LatticeMapping solution such that:
///
///     F * L1 * T * N = L2
///
/// Therefore, the AtomMapping translation is:
///
///     translation = F * trial_translation
///
///
/// \param assignment Assignment solution, with the convention
///     `atom_index = assignment[site_index]`. With the cost
///     matrix constructed according to the convention
///     `cost_matrix(site_index, atom_index)` this is equal
///     to the AtomMapping permutation.
/// \param site_displacements The site-to-atom displacements,
///     of minimum length under periodic boundary conditions
///     as used in the assignment problem.
/// \param trial_translation A translation applied to atom
///     coordinates to bring the atoms and sites into alignment.
/// \param deformation_gradient The deformation gradient of
///     the lattice mapping that this atom mapping is solved
///     in the context of.
///
/// \returns An AtomMapping constructed from the inputs, with
///     mean displacement removed.
AtomMapping make_atom_mapping_from_assignment(
    murty::Assignment const &assignment,
    std::vector<std::vector<Eigen::Vector3d>> const &site_displacements,
    Eigen::VectorXd trial_translation,
    Eigen::Matrix3d const &deformation_gradient) {
  Index N_site = assignment.size();

  // atom[perm[i]] -> is assigned to -> site[i]
  auto const &perm = assignment;

  // get displacements
  Eigen::MatrixXd disp(3, N_site);
  for (Index site_index = 0; site_index < N_site; ++site_index) {
    disp.col(site_index) = site_displacements[site_index][perm[site_index]];
  }

  // calculate mean displacement
  Eigen::Vector3d mean_disp = Eigen::Vector3d::Zero();
  for (Index site_index = 0; site_index < N_site; ++site_index) {
    mean_disp += disp.col(site_index);
  }
  mean_disp /= double(N_site);

  // subtract mean displacement
  for (Index site_index = 0; site_index < N_site; ++site_index) {
    disp.col(site_index) -= mean_disp;
  }
  trial_translation -= mean_disp;

  return AtomMapping(disp, perm, deformation_gradient * trial_translation);
}

}  // namespace mapping_impl

/// \brief Constructor
MappingNode::MappingNode(
    double _lattice_cost,
    std::shared_ptr<LatticeMappingSearchData const> _lattice_mapping_data,
    double _atom_cost,
    std::shared_ptr<AtomMappingSearchData const> _atom_mapping_data,
    murty::Node _assignment_node, AtomMapping _atom_mapping, double _total_cost)
    : lattice_cost(_lattice_cost),
      lattice_mapping_data(std::move(_lattice_mapping_data)),
      atom_cost(_atom_cost),
      atom_mapping_data(std::move(_atom_mapping_data)),
      assignment_node(std::move(_assignment_node)),
      atom_mapping(std::move(_atom_mapping)),
      total_cost(_total_cost) {}

/// \struct MappingSearch
/// \brief Performs structure mapping searches
///
/// The MappingSearch structure includes parameters, data,
/// and methods used to search for low cost structure mappings.
///
/// It holds a queue of MappingNode, which encode a particular
/// structure mapping, and the data necessary to start from
/// that structure mapping and find sub-optimal atom mappings
/// as part of a search using the Murty Algorithm for
/// sub-optimal assignments.
///
/// It also holds a sorted container of the best results found
/// so far which satisfy some acceptance criteria:
/// - min_cost: Keep mappings with total cost >= min_cost
/// - max_cost: Keep mappings with total cost <= max_cost
/// - k_best: Keep the k_best mappings with lowest total cost
///   that also satisfy the min/max cost criteria.
///
/// It also holds an `overflow` container to keep approximate
/// ties with the k_best result.
///
/// Overview of methods:
/// - `MappingSearch::make_and_insert_mapping_node`:
///   - Given a lattice mapping, trial translation, and optionally assignments
///   to force on or off, this solves the assignment problem, constructs a
///   MappingNode holding the solution. The MappingNode is inserted in
///   this->queue (always) and this->results (if it satisfies the acceptance
///   criteria).
///   - This is typically used to initialize a search from one or more lattice
///   mappings and trial translations
/// - `MappingSearch::partition`:
///   - Given a MappingNode, this applies to Murty Algorithm to find the next
///   level of sub-optimal assignments, and for each constructs a MappingNode
///   holding the solution, and inserts the MappingNode in this->queue (always)
///   and this->results (if it satisfies the acceptance criteria).
///   - This is typically used to continue a search from the first element of
///   the queue. After partitioning, the element can be erased from the queue.
///
/// For more details, see J.C. Thomas, A.R. Natarajan, and
/// A.V. Van der Ven, npj Computational Materials (2021)7:164;
/// https://doi.org/10.1038/s41524-021-00627-0
///

/// \brief Constructor
///
/// \param _min_cost Keep mappings with total cost >= min_cost
/// \param _max_cost Keep mappings with total cost <= max_cost
/// \param _lattice_cost_weight The fraction of the total cost due to
///     the lattice strain cost. The remaining fraction
///     (1.-lattice_cost_weight) is due to the atom cost. Default=0.5.
/// \param _atom_cost_method One of "isotropic_atom_cost" or
///     "symmetry_breaking_atom_cost"
/// \param _k_best Keep the k_best results satisfying the min_cost and
///     max_cost constraints. If there are approximate ties, those
///     will also be kept.
/// \param _infinity The value to use in the assignment problem cost
///     matrix for unallowed assignments
/// \param _cost_tol Tolerance for checking if mapping costs are
///     approximately equal
MappingSearch::MappingSearch(double _min_cost, double _max_cost,
                             double _lattice_cost_weight,
                             std::string _atom_cost_method, int _k_best,
                             double _infinity, double _cost_tol)
    : min_cost(_min_cost),
      max_cost(_max_cost),
      lattice_cost_weight(_lattice_cost_weight),
      is_symmetry_breaking_atom_cost(
          mapping_impl::make_is_symmetry_breaking_atom_cost(_atom_cost_method)),
      k_best(_k_best),
      infinity(_infinity),
      cost_tol(_cost_tol) {}

double MappingSearch::make_atom_cost(
    Eigen::MatrixXd const &displacement,
    LatticeMappingSearchData const &lattice_mapping_data) const {
  auto const &prim_data = *(lattice_mapping_data.prim_data);
  auto const &L1 = prim_data.prim_lattice.lat_column_mat();
  auto const &lattice_mapping = lattice_mapping_data.lattice_mapping;

  if (this->is_symmetry_breaking_atom_cost) {
    if (!prim_data.prim_sym_invariant_displacement_modes.has_value()) {
      throw std::runtime_error(
          "Error in MappingSearch::make_atom_cost: The prim symmetry invariant "
          "displacement modes do not exist");
    }
    return make_symmetry_breaking_atom_cost(
        L1, lattice_mapping, displacement,
        lattice_mapping_data.unitcellcoord_index_converter,
        *(prim_data.prim_sym_invariant_displacement_modes));
  } else {
    return make_isotropic_atom_cost(L1, lattice_mapping, displacement);
  }
}

double MappingSearch::make_total_cost(double lattice_cost,
                                      double atom_cost) const {
  return this->lattice_cost_weight * lattice_cost +
         (1. - this->lattice_cost_weight) * atom_cost;
};

/// \brief Make a new MappingNode from an assignment problem node
///     with a solved sub_assignment
///
/// When constructing the AtomMapping component of a MappingNode,
/// this removes mean displacements and makes the proper
/// AtomMapping displacements and translation. It calculates
/// the atom_cost and total_cost using the parameters specified
/// at MappingSearch construction time.
MappingNode MappingSearch::make_mapping_node_from_assignment_node(
    murty::Node assignment_node, double lattice_cost,
    std::shared_ptr<LatticeMappingSearchData const> lattice_mapping_data,
    std::shared_ptr<AtomMappingSearchData const> atom_mapping_data) {
  AtomMapping atom_mapping = mapping_impl::make_atom_mapping_from_assignment(
      murty::make_assignment(assignment_node),
      atom_mapping_data->site_displacements,
      atom_mapping_data->trial_translation_cart,
      lattice_mapping_data->lattice_mapping.deformation_gradient);
  double atom_cost =
      this->make_atom_cost(atom_mapping.displacement, *lattice_mapping_data);
  double total_cost = this->make_total_cost(lattice_cost, atom_cost);
  return MappingNode(lattice_cost, std::move(lattice_mapping_data), atom_cost,
                     std::move(atom_mapping_data), std::move(assignment_node),
                     std::move(atom_mapping), total_cost);
}

/// \brief Insert mapping node into this->queue & this->results,
///     maintaining k-best results
///
/// A MappingNode is inserted into this->queue (always!), and then
/// inserted into this->results if it satisfies the min/max cost
/// and k-best criteria. Approximate ties with the k-best cost
/// are kept in this->overflow.
///
/// \param mapping_node MappingNode to insert
///
/// \returns Iterator to MappingNode in queue
///
std::multiset<MappingNode>::iterator MappingSearch::insert(
    MappingNode mapping_node) {
  std::multiset<MappingNode>::iterator it =
      this->queue.insert(std::move(mapping_node));

  // --- maintain k_best results, keeping ties in overflow ---
  // note: max_cost is modified by `maintain_k_best_results` to shrink
  // to the current k_best-th cost once k_best results are found
  if (it->total_cost > (this->min_cost - this->cost_tol)) {
    if (results.size() < k_best ||
        it->total_cost < this->max_cost + this->cost_tol) {
      this->results.emplace(
          StructureMappingCost(it->lattice_cost, it->atom_cost, it->total_cost),
          StructureMapping(it->lattice_mapping_data->prim_data->shared_prim,
                           it->lattice_mapping_data->lattice_mapping,
                           it->atom_mapping));
      mapping::maintain_k_best_results(
          this->k_best, this->cost_tol, this->results, this->overflow,
          [](StructureMappingCost const &key) { return key.total_cost; });
    }
  }

  return it;
}

/// \brief Make assignment and insert mapping node
///     into this->queue & this->results, maintaining k-best results
///
/// The (constrained) assignment problem is solved in context of
/// a particular lattice mapping and trial translation, and
/// the resulting AtomMapping, atom mapping cost, and total cost
/// are stored in a MappingNode, along with the (constrained)
/// assignment_node which allows continuing the search for
/// suboptimal solutions. The MappingNode is inserted in
/// this->queue. It is also inserted in this->results, if it
/// satisifies the cost range and k-best criteria.
///
///
/// \param lattice_cost The lattice mapping cost
/// \param lattice_mapping_data Data associated with the lattice
///     mapping that is the context in which the atom mapping is done
/// \param atom_mapping_data Data that can be used for all
///     atom mappings with the same lattice mapping and trial translation
/// \param forced_on A map of {site_index, atom_index} of
///     assignments that are forced on
/// \param forced_off A vector of {site_index, atom_index} of
///     assignments that are forced off (given infinity cost)
///
/// \returns Iterator to resulting MappingNode in queue
std::multiset<MappingNode>::iterator
MappingSearch::make_and_insert_mapping_node(
    double lattice_cost,
    std::shared_ptr<LatticeMappingSearchData const> lattice_mapping_data,
    std::shared_ptr<AtomMappingSearchData const> atom_mapping_data,
    std::map<Index, Index> forced_on,
    std::vector<std::pair<Index, Index>> forced_off) {
  // --- Find optimal assignment ---
  murty::Node assignment_node =
      murty::make_node(atom_mapping_data->cost_matrix, std::move(forced_on),
                       std::move(forced_off));
  std::tie(assignment_node.cost, assignment_node.sub_assignment) =
      murty::make_sub_assignment(
          hungarian::solve, atom_mapping_data->cost_matrix,
          assignment_node.unassigned_rows, assignment_node.unassigned_cols,
          assignment_node.forced_off, this->infinity, this->cost_tol);

  // --- Make mapping node from assignment solution ---
  MappingNode mapping_node = this->make_mapping_node_from_assignment_node(
      std::move(assignment_node), lattice_cost, std::move(lattice_mapping_data),
      atom_mapping_data);

  // --- Insert mapping node in queue and results, return queue iterator ---
  return this->insert(std::move(mapping_node));
}

/// \brief Make the next level of sub-optimal assignments and
///     inserts them into this->queue & this->results, maintaining
///     k-best results
///
/// The Murty algorithm is used to generate sub-optimal assignments
/// from a previous assignment solution. The resulting MappingNode
/// are inserted in this->queue. They are also inserted in
/// this->results, if they satisify the cost range and k-best
/// criteria.
///
/// \param node A MappingNode to find sub-optimal assignments from
///     using the Murty algorithm
///
/// \returns A vector of iterators to the generated sub-nodes in queue
///
std::vector<std::multiset<MappingNode>::iterator> MappingSearch::partition(
    MappingNode const &node) {
  // -- Make the next level of sub-optimal assignment solutions ---
  std::multiset<murty::Node> s;
  murty::partition(s, hungarian::solve, node.atom_mapping_data->cost_matrix,
                   node.assignment_node, this->infinity, this->cost_tol);

  // The sub-optimal assignment solutions are in 's',
  // and we want them to all end up in MappingNode.
  // This loop extracts the values from 's', uses them to
  // construct MappingNode, and continues until they are
  // all extracted.
  std::vector<std::multiset<MappingNode>::iterator> result;
  while (s.size()) {
    // --- Make mapping node from sub-optimal assignment ---
    MappingNode mapping_node = this->make_mapping_node_from_assignment_node(
        std::move(s.extract(s.begin()).value()), node.lattice_cost,
        node.lattice_mapping_data, node.atom_mapping_data);

    // --- Insert mapping node in queue and results, return queue iterator ---
    result.emplace_back(this->insert(std::move(mapping_node)));
  }

  return result;
}

}  // namespace mapping
}  // namespace CASM
