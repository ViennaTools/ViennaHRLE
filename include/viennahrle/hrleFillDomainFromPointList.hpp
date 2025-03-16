#ifndef HRLE_DOMAIN_FROM_POINT_LIST_HPP
#define HRLE_DOMAIN_FROM_POINT_LIST_HPP

#include <algorithm>
#include <vector>

#include "hrleDomain.hpp"

namespace viennahrle {
using namespace viennacore;

// this functions clears the hrleDomain
// and fills it with the Data in pointData,
// a sorted list of index/value pairs
template <class T, int D>
void FillDomainFromPointList(
    Domain<T, D> &newDomain,
    std::vector<std::pair<VectorType<IndexType, D>, T>> pointData,
    const T &backgroundValue, const bool sortPointList = true) {

  typedef std::pair<VectorType<IndexType, D>, T> indexValuePairType;

  // TODO: is not parallelized yet
  newDomain.initialize();

  // if there are no points, just initialize an empty hrleDomain
  if (pointData.empty()) {
    return;
  }

  if (sortPointList) {
    std::sort(
        pointData.begin(), pointData.end(),
        [](const indexValuePairType &a, const indexValuePairType &b) -> bool {
          return a.first < b.first;
        });
  }

  const auto &grid = newDomain.getGrid();

  if (pointData.front().first != grid.getMinGridPoint()) {
    newDomain.insertNextUndefinedPoint(0, grid.getMinGridPoint(),
                                       backgroundValue);
  }

  auto pointDataBegin = pointData.begin(); //+starts[t_num];
  auto pointDataEnd = pointData.end();     // pointData.begin()+starts[t_num+1];

  auto pointDataIt = pointDataBegin;
  while (pointDataIt != pointDataEnd) {

    // Add defined point as it appears in the list
    newDomain.insertNextDefinedPoint(0, pointDataIt->first,
                                     pointDataIt->second);

    VectorType<IndexType, D> index = pointDataIt->first;
    VectorType<IndexType, D> next_index;

    ++pointDataIt;

    // choose correct next index
    if (pointDataIt == pointDataEnd) {
      next_index = grid.getMaxGridPoint();
      ++next_index[D - 1];
    } else {
      next_index = pointDataIt->first;
    }

    for (int q = 0; q < D; q++) {
      VectorType<IndexType, D> tmp = index;
      ++tmp[q];
      if (tmp[q] > grid.getMaxGridPoint(q))
        continue;
      for (int r = 0; r < q; ++r)
        tmp[r] = grid.getMinGridPoint(r);

      if (tmp >= next_index)
        break;

      newDomain.insertNextUndefinedPoint(0, tmp, backgroundValue);
    }
  }

  newDomain.finalize();
}
} // namespace viennahrle

#endif // HRLE_DOMAIN_FROM_POINT_LIST_HPP
