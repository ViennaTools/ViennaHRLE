#ifndef HRLE_FOMAIN_FROM_DATA_HPP
#define HRLE_FOMAIN_FROM_DATA_HPP

#include <algorithm>
#include <vector>

#include "hrleDomain.hpp"

// this functions clears the hrleDomain
// and fills it with the Data in pointData,
// a sorted list of index/value pairs
template <class T, int D>
void hrleFillDomainFromPointList(
    hrleDomain<T, D> &newDomain,
    std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>> pointData,
    const T &backgroundValue, const bool sortPointList = true) {

  typedef std::pair<hrleVectorType<hrleIndexType, D>, T> indexValuePairType;

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

  const hrleGrid<D> &grid = newDomain.getGrid();

  if (pointData.front().first != grid.getMinIndex()) {
    newDomain.insertNextUndefinedPoint(0, grid.getMinIndex(), backgroundValue);
  }

  auto pointDataBegin = pointData.begin(); //+starts[t_num];
  auto pointDataEnd = pointData.end();     // pointData.begin()+starts[t_num+1];

  hrleVectorType<hrleIndexType, D> currentIndex = pointDataBegin->first;

  auto pointDataIt = pointDataBegin;
  while (pointDataIt != pointDataEnd) {

    // Add defined point as it appears in the list
    newDomain.insertNextDefinedPoint(0, pointDataIt->first,
                                     pointDataIt->second);

    hrleVectorType<hrleIndexType, D> index = pointDataIt->first;
    hrleVectorType<hrleIndexType, D> next_index;

    pointDataIt++;

    // choose correct next index
    if (pointDataIt == pointDataEnd) {
      next_index = grid.getMaxIndex();
      next_index[D - 1]++;
    } else {
      next_index = pointDataIt->first;
    }

    for (int q = 0; q < D; q++) {
      hrleVectorType<hrleIndexType, D> tmp = index;
      tmp[q]++;
      if (tmp[q] > grid.getMaxIndex(q))
        continue;
      for (int r = 0; r < q; ++r)
        tmp[r] = grid.getMinIndex(r);

      if (tmp >= next_index)
        break;

      newDomain.insertNextUndefinedPoint(0, tmp, backgroundValue);
    }
  }

  newDomain.finalize(2);
}

#endif // HRLE_FOMAIN_FROM_DATA_HPP
