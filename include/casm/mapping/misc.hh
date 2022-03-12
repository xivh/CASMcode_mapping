#ifndef CASM_mapping_misc
#define CASM_mapping_misc

#include <map>
#include <optional>

namespace CASM {
namespace mapping {

/// \brief Maintain k-best results
///
/// \param k_best The optional number of solutions, T, to keep. Solutions
///     approximately equal to the k-th best costing solution are also
///     kept, in `overflow`.
/// \param cost_tol Tolerance used for comparing costs
/// \param results A map of cost -> solution
/// \param overflow A map of cost -> solution for solutions approximately equal
///     to results.rbegin()
/// \param max_cost The max cost of solution to keep in results (kept if cost
///     < max_cost + cost_tol). If results.size() > *k_best, the max_cost is
///     shrunk to results.rbegin()->first.
///
/// This does the following:
/// - Nothing if !k_best.has_value()
/// - While results.size() > *k_best:
///   - Move results approximately equal to results.rbegin() to overflow, until
///     results.size() == *k_best
///   - Shrink max_cost to equal results.rbegin()->first
template <typename K, typename T>
void maintain_k_best_results(std::optional<int> const &k_best, double cost_tol,
                             std::multimap<K, T> &results,
                             std::multimap<K, T> &overflow, double &max_cost) {
  if (!k_best.has_value()) {
    return;
  }
  // only keep k-best results, but include those ~equal to k-th best result
  while (static_cast<int>(results.size()) > *k_best) {
    // check if last two results are approximately equal
    // reverse iterator, points at extra node
    auto last_it = --results.end();

    // reverse iterator, points at k-th best result
    auto nexttolast_it = --(--results.end());

    // if extra is not approximately equal cost to k-th best, erase extra and
    // overflow
    if ((last_it->first.cost - nexttolast_it->first.cost) > cost_tol) {
      results.erase(last_it);
      overflow.clear();
    }
    // if extra is approximately equal to last, move extra to overflow
    else {
      overflow.insert(results.extract(last_it));
    }
    // shrink max_cost
    max_cost = results.rbegin()->first.cost;
  }
}

}  // namespace mapping
}  // namespace CASM

#endif
