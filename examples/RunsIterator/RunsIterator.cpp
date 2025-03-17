#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  using namespace viennahrle;

  // set domain bounds
  IndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -20;
    max[i] = 20;
  }

  // declare domain with bounds min,max
  Grid<D> grid(min, max);
  Domain<char, D> data(grid);

  std::vector<std::pair<Index<D>, char>> pointData;
  Index<D> index(0, 5, 0);

  std::string helloString = "Hello, World!";

  for (char &i : helloString) {
    pointData.emplace_back(index, i);
    index[0] += 2;
    index[1] += 1;
  }

  FillDomainFromPointList(data, pointData,
                          '.'); // last parameter is the background value to use

  // iterate over hrle structure and output the values
  SparseIterator<Domain<char, D>> it(data);
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
