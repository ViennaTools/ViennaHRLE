#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseIterator.hpp>
#include <hrleSparseOffsetIterator.hpp>
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

  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 1;
    // index[1] += 1;
  }
  ++index[1];
  index[0] = 0;
  for (int i = int(helloString.size()) - 1; i >= 0; --i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 1;
    // index[1] += 1;
  }

  FillDomainFromPointList(data, pointData,
                          '.'); // last parameter is the background value to use

  // iterate over  structure and output each value and the one above it
  Index<D> offset(0, -1, 0);
  SparseIterator<Domain<char, D>> it(data);
  SparseOffsetIterator<Domain<char, D>> offsetIt(data, offset);
  while (!it.isFinished()) {
    std::cout << it.getStartIndices() << " " << it.getValue() << "\t\t"
              << offsetIt.getStartIndices() - offset << " "
              << offsetIt.getValue() << std::endl;
    it.next();
    if (Compare(it.getStartIndices(), offsetIt.getStartIndices()) > 0)
      offsetIt.next();
  }

  // print data structure in debug output
  // data.print();

  return 0;
}
