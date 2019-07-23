#include <hrleFillDomainFromPointList.hpp>
#include <hrleIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;

  // set domain bounds
  hrleIndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -5;
    max[i] = 25;
  }

  // declare domain with bounds min,max
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index(0);

  std::string helloString = "Hello, World!";

  for (hrleIndexType i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  hrleFillDomainFromPointList(
      data, pointData,
      '.'); // last parameter is the background value to use

  // iterate over hrle structure and output the values showing the dense data
  // set
  std::cout << "Dense Data Set filled with background value:" << std::endl;
  hrleConstIterator<hrleDomain<char, D>> it(data);
  int y = data.getGrid().getMinIndex(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      std::cout << std::endl;
    }
    std::cout << (it++).getValue() << " ";
  }

  std::cout << it.getValue() << std::endl;

  return 0;
}
