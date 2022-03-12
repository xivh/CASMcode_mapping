#ifndef CASM_mapping_murty
#define CASM_mapping_murty

#include <functional>
#include <vector>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace mapping {
namespace murty {

/// \brief The assignment vector gives `j = assignment[i]`, where i is the
///     "worker" (row) and j is the assigned "task" (column).
typedef std::vector<Index> Assignment;

/// \brief The method used to solve for an optimal assignment
///
/// Expected parameters:
/// - Eigen::MatrixXd const &cost_matrix
/// - double infinity
/// - double tol
///
/// Expected results:
/// - {cost, assignment}, if optimal solution is found
/// - {infinity, {}}, if could not solve (no solution that doesn't
///   include a forced_off assignment)
typedef std::function<std::pair<double, Assignment>(Eigen::MatrixXd const &,
                                                    double, double)>
    AssignmentMethod;

/// \brief Find the k best solutions to the assignment problem
///     using the Murty algorithm
std::vector<std::pair<double, Assignment>> solve(
    AssignmentMethod assign_f, Eigen::MatrixXd const &cost_matrix, int k_best,
    double min_cost = 0.0, double max_cost = 1e20, double infinity = 1e20,
    double tol = 1e-5);

}  // namespace murty
}  // namespace mapping
}  // namespace CASM

#endif
