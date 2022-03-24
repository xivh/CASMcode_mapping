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
///
/// This does the following:
/// - Nothing if !k_best.has_value()
/// - While results.size() > *k_best:
///   - Move results approximately equal to results.rbegin() to overflow, until
///     results.size() == *k_best
template <typename K, typename T, typename Compare, typename GetCostFromKey>
void maintain_k_best_results(std::optional<int> const &k_best, double cost_tol,
                             std::multimap<K, T, Compare> &results,
                             std::multimap<K, T, Compare> &overflow,
                             GetCostFromKey get_cost_f) {
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
    if ((get_cost_f(last_it->first) - get_cost_f(nexttolast_it->first)) >
        cost_tol) {
      results.erase(last_it);
      overflow.clear();
    }
    // if extra is approximately equal to last, move extra to overflow
    else {
      overflow.insert(results.extract(last_it));
    }
  }
}

}  // namespace mapping
}  // namespace CASM

#endif
