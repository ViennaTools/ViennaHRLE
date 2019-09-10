#ifndef HRLE_FILL_DOMAIN_WITH_SIGNED_DISTANCE_HPP
#define HRLE_FILL_DOMAIN_WITH_SIGNED_DISTANCE_HPP

/// Helper function to fill hrleDomain with narrowband/sparse-field
/// signed distance function values. Empty space (undefined runs) is
/// filled with the sign of the last defined run.
//
/// NOTE: a valid list of points must include all grid points, which are
/// connected by edges of the grid, which are intersected by the surface
template <class T, int D>
void hrleFillDomainWithSignedDistance(
    hrleDomain<T, D> &newDomain,
    std::vector<std::pair<hrleVectorType<hrleIndexType, D>, T>> pointData,
    const T &negValue, const T &posValue, const bool sortPointList = true) {

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
    newDomain.insertNextUndefinedPoint(
        0, grid.getMinIndex(),
        (pointData.front().second < 0) ? negValue : posValue);
  }

  auto pointDataBegin = pointData.begin(); //+starts[t_num];
  auto pointDataEnd = pointData.end();     // pointData.begin()+starts[t_num+1];

  hrleVectorType<hrleIndexType, D> currentIndex = pointDataBegin->first;

  auto pointDataIt = pointDataBegin;

  hrleVectorType<bool, D> signs(pointDataBegin->second < 0);

  while (pointDataIt != pointDataEnd) {
    bool setPoint = true;

    // if boundary conditions are infinite always set the point
    // if not, check, whether it is inside of domain
    for (unsigned i = 0; i < D; ++i) {
      if (grid.getBoundaryConditions(i) != hrleGrid<D>::INFINITE_BOUNDARY) {
        if (pointDataIt->first[i] > grid.getMaxBounds(i) ||
            pointDataIt->first < grid.getMinBounds())
          setPoint = false;
      }
    }

    if (setPoint) {
      // Add defined point as it appears in the list
      newDomain.insertNextDefinedPoint(0, pointDataIt->first,
                                       pointDataIt->second);

      // determine signs for next undefined runs
      {
        bool changeSign = false;
        for (int i = D - 1; i >= 0; --i) {
          changeSign = changeSign || (pointDataIt->first[i] > currentIndex[i]);
          if (changeSign) {
            signs[i] = pointDataIt->second < 0;
            currentIndex[i] = pointDataIt->first[i];
          }
        }
      }
    }

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

    // move current index by one grid spacing and see if the next
    // point has the same index, if not, there must be an undefined
    // run inbetween
    for (int q = 0; q < D; q++) {
      hrleVectorType<hrleIndexType, D> tmp = index;
      tmp[q]++;
      if (tmp[q] > grid.getMaxIndex(q))
        continue;
      for (int r = 0; r < q; ++r)
        tmp[r] = grid.getMinIndex(r);

      if (tmp >= next_index)
        break;

      newDomain.insertNextUndefinedPoint(0, tmp,
                                         signs[q] ? negValue : posValue);
    }
  }

  newDomain.finalize();
}

#endif // HRLE_FILL_DOMAIN_WITH_SIGNED_DISTANCE_HPP
