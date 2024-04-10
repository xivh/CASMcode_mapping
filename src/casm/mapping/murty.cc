#include "casm/mapping/murty.hh"

namespace CASM {
namespace mapping {
namespace murty {
namespace murty_impl {

/// \brief Returns the cost due to just "forced_on" assignments
double make_forced_on_cost(Eigen::MatrixXd const &cost_matrix,
                           Node const &node) {
  double forced_on_cost = 0.0;
  for (auto const &pair : node.forced_on) {
    forced_on_cost += cost_matrix(pair.first, pair.second);
  }
  return forced_on_cost;
}

}  // namespace murty_impl

/// \brief Find the k best solutions to the assignment problem
///     using the Murty algorithm
///
/// The assignment problem is: minimize the cost of assigning m
/// "workers" to n "tasks", where the cost of assigning "worker" i
/// to "task" j is cost_matrix(i,j).
///
/// The Murty algorithm finds the optimal solution, then checks
/// the solutions of constrained sub-problems where some of the
/// optimal assignments are forced on (the assignment is made
/// in all solutions considered) and some of the optimal
/// assignments are forced "off" (the cost is set to infinity
/// so assignment is never made). This is done in a systematic
/// fashion forming a "best-first search" that can find all
/// solutions.
///
/// This implementation allows stopping the search at a certain
/// number (k) of solutions, or at a certain maximum cost. The
/// solutions k+1, k+2, etc. which are equivalent to the k-best
/// solution or maximum cost within a specified tolerance will
/// also be returned.
///
/// \param assign_f Method used for calculating the best assignment
/// \param cost_matrix The cost of assigning "worker" i to "task" j
///     is cost_matrix(i,j). The number of rows and columns must
///     be greater than 1. The number of rows must be equal to the
///     number of columns.
/// \param k_best Number of solutions to return. Must be greater
///     than or equal to 1.
/// \param min_cost Minimum cost solutions to include in the results.
///     If nullopt, the min_cost value is set to the minimum
///     coefficient value.
/// \param max_cost Maximum cost solutions to include in the results.
///     If nullopt, the max_cost value is set to the infinity value.
/// \param infinity Cost used for "infinity", when an assignment is
///     forced off
/// \param tol Tolerance used for comparing costs
///
/// \returns A vector of best solutions, as pairs of {cost, assignment}.
///     The first solution is the optimal cost and assignment. The
///     assignment vector gives `j = assignment[i]`, where i is the
///     worker (row) and j is the task (column). The cost of all
///     solutions will be less than max_cost + tol. More than k
///     solutions may be returned if solutions k+1, k+2, etc. have
///     cost less than the cost of the k-th solution + tol.
std::vector<std::pair<double, Assignment>> solve(
    AssignmentMethod assign_f, Eigen::MatrixXd const &cost_matrix, int k_best,
    std::optional<double> min_cost, std::optional<double> max_cost,
    double infinity, double tol) {
  // --- Input validation ---
  if (k_best < 1) {
    throw std::runtime_error("Error in murty::solve: k_best < 1");
  }
  if (cost_matrix.rows() < 1) {
    throw std::runtime_error("Error in murty::solve: cost_matrix.rows() < 1");
  }
  if (cost_matrix.cols() < 1) {
    throw std::runtime_error("Error in murty::solve: cost_matrix.cols() < 1");
  }
  if (cost_matrix.rows() != cost_matrix.cols()) {
    throw std::runtime_error(
        "Error in murty::solve: cost_matrix.rows() != cost_matrix.cols()");
  }

  // --- Defaults ---
  if (!min_cost.has_value()) {
    min_cost = cost_matrix.minCoeff();
  }
  if (!max_cost.has_value()) {
    max_cost = infinity;
  }

  std::vector<std::pair<double, Assignment>> results;

  // --- Find optimal assignment ---
  Node optimal_node = make_node(cost_matrix);
  std::tie(optimal_node.cost, optimal_node.sub_assignment) =
      make_sub_assignment(assign_f, cost_matrix, optimal_node.unassigned_rows,
                          optimal_node.unassigned_cols, optimal_node.forced_off,
                          infinity, tol);

  // If cost is infinite or cost is greater than max_cost; then return empty
  // results
  if ((infinity - tol < optimal_node.cost) ||
      (*max_cost + tol <= optimal_node.cost)) {
    return results;
  }
  // If cost is greater than or equal to min_cost; then add it to results
  if (*min_cost - tol < optimal_node.cost) {
    results.emplace_back(optimal_node.cost, make_assignment(optimal_node));
  }

  // --- Find suboptimal assignments ---
  std::multiset<Node> nodes;
  nodes.insert(optimal_node);

  // Get iterator to current best node
  auto node_it = nodes.begin();
  while (true) {
    // Partition by current best node
    partition(nodes, assign_f, cost_matrix, *node_it, infinity, tol);

    // Get next best node and its cost
    nodes.erase(node_it);
    if (nodes.size() == 0) {
      break;
    }
    node_it = nodes.begin();

    // If cost is less than min_cost, continue
    double cost = node_it->cost;
    if (cost <= *min_cost - tol) {
      continue;
    }
    // If cost is infinite or cost is greater than max_cost; then return results
    if ((infinity - tol < cost) || (*max_cost + tol <= cost)) {
      break;
    }
    // If k_best not satisfied or tied with k_best result; then add to results
    if (results.size() < k_best || cost < results[k_best - 1].first + tol) {
      results.emplace_back(cost, make_assignment(*node_it));
      continue;
    }
    // Otherwise, k_best is satisfied and the current cost is not tied
    // with the k_best result; so break
    break;
  }
  return results;
}

/// \brief Returns a Node representing the (constrained) assignment problem
///
/// \param cost_matrix The cost matrix for the assignement problem
/// \param forced_on A map of indicies of {row, column} of
///     assignments that are forced on
/// \param forced_off A vector of indicies of {row, column} of
///     assignments that are forced off (given infinity cost)
///
/// \returns node The resulting Node has default initialized sub_assignment
///     and cost, while unassigned_rows and unassigned_cols are set to be
///     consistent with forced_on.
Node make_node(Eigen::MatrixXd const &cost_matrix,
               std::map<Index, Index> forced_on,
               std::vector<std::pair<Index, Index>> forced_off) {
  Node node;
  for (Index i = 0; i < cost_matrix.rows(); ++i) {
    node.unassigned_rows.insert(i);
  }
  for (Index i = 0; i < cost_matrix.cols(); ++i) {
    node.unassigned_cols.insert(i);
  }
  for (auto const &pair : forced_on) {
    node.unassigned_rows.erase(pair.first);
    node.unassigned_cols.erase(pair.second);
  }
  node.forced_on = std::move(forced_on);
  node.forced_off = std::move(forced_off);

  return node;
}

/// \brief Returns the full assignment
///
/// This combines forced_on and sub_assignment and copies
/// the solution to a vector, under the assumption that
/// all workers must be assigned.
Assignment make_assignment(Node const &node) {
  std::map<Index, Index> ordered = node.forced_on;
  for (auto const &x : node.sub_assignment) {
    ordered.emplace(x);
  }
  Assignment assignment;
  for (auto const &x : ordered) {
    assignment.push_back(x.second);
  }
  return assignment;
}

/// \brief Return the cost for a particular assignment
double make_cost(Eigen::MatrixXd const &cost_matrix,
                 Assignment const &assignment) {
  double cost = 0.0;
  for (Index i = 0; i < assignment.size(); ++i) {
    cost += cost_matrix(i, assignment[i]);
  }
  return cost;
}

/// \brief Partition a Node, adding results to a multiset of Node
///
/// \param node_set A multiset of Node, sorted by ascending cost, containing
///     suboptimal assignments, as constructed by Murty Algorithm partitioning
/// \param assign_f Method used for calculating the best assignment
/// \param cost_matrix The cost of assigning "worker" i to "task" j is
///     cost_matrix(i,j)
/// \param node Encodes the "current best" assignment
/// \param infinity Cost used for "infinity", when an assignment is forced
///     off
/// \param tol Tolerance used for comparing costs
///
/// \returns A Node encoding the "next best" assignment
void partition(std::multiset<Node> &node_set, AssignmentMethod assign_f,
               Eigen::MatrixXd const &cost_matrix, Node const &node,
               double infinity, double tol) {
  using namespace murty_impl;
  // node partitioning:
  // - where node.sub_assignment == {x0, x1, x2, x3, ...}, xi={row,column}
  // sub-node 0: {node.forced_on, node.forced_off + x0}
  // sub-node 1: {node.forced_on + x0, node.forced_off + x1}
  // sub-node 2: {node.forced_on + x0 + x1, node.forced_off + x2}
  // sub-node 3: {node.forced_on + x0 + x1 + x2, node.forced_off + x3}
  // ...

  if (node.unassigned_rows.size() == 1) {
    // no sub-assignments in this case
    return;
  }

  Node subnode = node;

  // 'x' is a particular assignement {worker/row, task/column}
  for (auto const &x : node.sub_assignment) {
    subnode.forced_off.push_back(x);

    double sub_assignment_cost;
    std::tie(sub_assignment_cost, subnode.sub_assignment) = make_sub_assignment(
        assign_f, cost_matrix, subnode.unassigned_rows, subnode.unassigned_cols,
        subnode.forced_off, infinity, tol);

    // handle failure to assign all "workers" (rows) by not saving the node
    if (subnode.sub_assignment.size() + subnode.forced_on.size() ==
        cost_matrix.rows()) {
      subnode.cost =
          sub_assignment_cost + make_forced_on_cost(cost_matrix, subnode);
      node_set.insert(subnode);
    }

    subnode.forced_on.emplace(x);
    if (subnode.forced_on.size() == cost_matrix.rows()) {
      break;
    }
    subnode.unassigned_rows.erase(x.first);
    subnode.unassigned_cols.erase(x.second);
    subnode.forced_off.pop_back();
  }
}

/// \brief Solve the assignment problem given certain assignments
///    forced on and certain assignments forced off
///
/// \param assign_f Method used for calculating the best assignment
/// \param cost_matrix The cost of assigning "worker" i to "task" j is
///     cost_matrix(i,j). This is the original, full, cost_matrix.
/// \param unassigned_rows Indices of "workers" (rows) not currently
///     forced into a particular assignment.
/// \param unassigned_rows Indices of "tasks" (columns) not currently
///     forced into a particular assignment.
/// \param forced_off Pairs of "workers" and "tasks" {row, column}
///     which are forced to be not be part of the assignment solution.
/// \param infinity Cost used for "infinity", when an assignment is forced
///     off
/// \param tol Tolerance used for comparing costs
///
/// \returns {sub_assignment_cost, sub_assignment}, where
///     sub_assignment_cost is the cost of the sub_assignment only,
///     and sub_assignment is the optimal assignment {row, column}
///      of the sub cost matrix formed from the unassigned_rows and
///     unassigned_cols of the cost matrix, given the constraints.
///
std::pair<double, std::map<Index, Index>> make_sub_assignment(
    AssignmentMethod assign_f, Eigen::MatrixXd cost_matrix,
    std::set<Index> const &unassigned_rows,
    std::set<Index> const &unassigned_cols,
    std::vector<std::pair<Index, Index>> const &forced_off, double infinity,
    double tol) {
  // --- Setup sub-assignment problem ---
  Eigen::MatrixXd sub_cost_matrix(unassigned_rows.size(),
                                  unassigned_cols.size());

  // Use infinity cost for assignments that are forced off
  for (auto const &pair : forced_off) {
    cost_matrix(pair.first, pair.second) = infinity;
  }

  // Use original cost for assignments that are not forced off
  Index row = 0;
  for (Index full_matrix_row : unassigned_rows) {
    Index col = 0;
    for (Index full_matrix_col : unassigned_cols) {
      sub_cost_matrix(row, col) = cost_matrix(full_matrix_row, full_matrix_col);
      ++col;
    }
    ++row;
  }

  // --- Solve sub-assignment problem ---
  std::pair<double, Assignment> tmp = assign_f(sub_cost_matrix, infinity, tol);
  std::pair<double, std::map<Index, Index>> result;
  result.first = tmp.first;

  // If cost is >= infinity...
  // no solution that doesn't include forced_off assignment
  if (result.first >= infinity) {
    return result;
  }

  // --- Convert sub-problem indices to full-problem indices ---
  auto unassigned_rows_it = unassigned_rows.begin();
  std::vector<Index> unassigned_cols_vec{unassigned_cols.begin(),
                                         unassigned_cols.end()};
  for (Index n : tmp.second) {
    result.second.emplace(*unassigned_rows_it, unassigned_cols_vec[n]);
    ++unassigned_rows_it;
  }
  return result;
}

}  // namespace murty
}  // namespace mapping
}  // namespace CASM
