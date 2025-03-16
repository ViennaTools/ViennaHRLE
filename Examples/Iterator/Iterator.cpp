#include <hrleDenseCellIterator.hpp>
#include <hrleDenseIterator.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseCellIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  using namespace viennahrle;

  // set domain bounds
  IndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -5;
    max[i] = 25;
  }

  // declare domain with bounds min,max
  Grid<D> grid(min, max);
  Domain<char, D> data(grid);

  std::vector<std::pair<Index<D>, char>> pointData;
  Index<D> index(0);

  std::string helloString = "Hello, World!";

  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  FillDomainFromPointList(data, pointData,
                          '.'); // last parameter is the background value to use

  // iterate over hrle structure and output the values showing the dense data
  // set
  std::cout << "Dense Data Set filled with background value:" << std::endl;
  ConstDenseIterator<Domain<char, D>> it(data);
  int y = data.getGrid().getMinIndex(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      std::cout << std::endl;
    }
    std::cout << (it++).getValue() << " ";
  }

  std::cout << std::endl;

  ConstDenseCellIterator<Domain<char, D>> it2(data);

  SparseCellIterator<Domain<char, D>> it3(data);

  return 0;
}
