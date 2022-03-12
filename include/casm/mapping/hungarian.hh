#ifndef CASM_mapping_hungarian
#define CASM_mapping_hungarian

#include <vector>

#include "casm/global/definitions.hh"
#include "casm/global/eigen.hh"

namespace CASM {
namespace mapping {
namespace hungarian {

/// \brief The assignment vector gives `j = assignment[i]`, where i is the
///     "worker" (row) and j is the assigned "task" (column).
typedef std::vector<Index> Assignment;

/// \brief Find the optimal solution to the assignment problem
std::pair<double, Assignment> solve(Eigen::MatrixXd const &cost_matrix,
                                    double infinity = 1e20, double tol = 1e-5);

}  // namespace hungarian
}  // namespace mapping
}  // namespace CASM

#endif
