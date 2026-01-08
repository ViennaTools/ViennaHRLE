#include <hrleDenseIterator.hpp>
#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <vcTestAsserts.hpp>
#include <vector>

using namespace viennahrle;

int main() {
  constexpr int D = 2;
  using DataType = int;

  IndexType minArr[D] = {-10, -10};
  IndexType maxArr[D] = {10, 10};

  Grid<D> grid(minArr, maxArr);
  Domain<DataType, D> domain(&grid);

  std::vector<std::pair<Index<D>, DataType>> points;
  points.emplace_back(Index<D>(0, 0), 10);
  points.emplace_back(Index<D>(1, 1), 20);
  points.emplace_back(Index<D>(-5, 5), 30);

  FillDomainFromPointList(domain, points, 0);

  ConstDenseIterator<Domain<DataType, D>> it(domain);

  int countDefined = 0;
  while (!it.isFinished()) {
    auto idx = it.getIndices();
    auto val = it.getValue();

    if (idx == Index<D>(0, 0)) {
      VC_TEST_ASSERT(val == 10);
      countDefined++;
    } else if (idx == Index<D>(1, 1)) {
      VC_TEST_ASSERT(val == 20);
      countDefined++;
    } else if (idx == Index<D>(-5, 5)) {
      VC_TEST_ASSERT(val == 30);
      countDefined++;
    } else {
      VC_TEST_ASSERT(val == 0); // Background value
    }
    ++it;
  }

  VC_TEST_ASSERT(countDefined == 3);
  std::cout << "Iterator consistency check passed." << std::endl;

  return 0;
}
