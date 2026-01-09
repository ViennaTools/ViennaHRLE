#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <vcTestAsserts.hpp>
#include <vector>

using namespace viennahrle;

int main() {
  constexpr int D = 3;
  using DataType = int;

  IndexType minArr[D];
  IndexType maxArr[D];
  for (int i = 0; i < D; ++i) {
    minArr[i] = -10;
    maxArr[i] = 10;
  }

  Grid<D> grid(minArr, maxArr);

  // Case 1: Solid block - contiguous defined points
  // Should result in fewer runs (geometry compression)
  SizeType runsSolid = 0;
  {
    Domain<DataType, D> domain(&grid);
    std::vector<std::pair<Index<D>, DataType>> points;

    // Fill a 5x5x5 block
    for (int z = -2; z <= 2; ++z) {
      for (int y = -2; y <= 2; ++y) {
        for (int x = -2; x <= 2; ++x) {
          points.emplace_back(Index<D>(x, y, z), 1);
        }
      }
    }

    FillDomainFromPointList(domain, points, 0);
    runsSolid = domain.getDomainSegment(0).getNumberOfRuns(0);
  }

  // Case 2: Checkerboard block - alternating defined/undefined
  // Should result in many more runs because of fragmentation
  SizeType runsSparse = 0;
  {
    Domain<DataType, D> domain(&grid);
    std::vector<std::pair<Index<D>, DataType>> points;

    // Fill a 5x5x5 block but skip every other point
    for (int z = -2; z <= 2; ++z) {
      for (int y = -2; y <= 2; ++y) {
        for (int x = -2; x <= 2; ++x) {
          if ((x + y + z) % 2 == 0) {
            points.emplace_back(Index<D>(x, y, z), 1);
          }
        }
      }
    }

    FillDomainFromPointList(domain, points, 0);
    runsSparse = domain.getDomainSegment(0).getNumberOfRuns(0);
  }

  std::cout << "Solid Runs: " << runsSolid
            << ", Sparse (Checkerboard) Runs: " << runsSparse << std::endl;

  VC_TEST_ASSERT(runsSolid < runsSparse);

  return 0;
}
