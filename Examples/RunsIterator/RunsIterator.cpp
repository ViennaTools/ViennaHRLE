#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;

  // set domain bounds
  hrleIndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -20;
    max[i] = 20;
  }

  // declare domain with bounds min,max
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index(0, 5, 0);

  std::string helloString = "Hello, World!";

  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  hrleFillDomainFromPointList(
      data, pointData,
      '.'); // last parameter is the background value to use

  // iterate over hrle structure and output the values
  hrleSparseIterator<hrleDomain<char, D>> it(data);
  while (!it.isFinished()) {
    // while iterating, change each 'o' to an 'a'
    if (it.getValue() == 'o')
      it.getValue() = 'a';

    std::cout << it.getStartIndices() << ": " << std::setw(10) << it.getValue()
              << std::endl;
    it.next();
  }

  return 0;
}
