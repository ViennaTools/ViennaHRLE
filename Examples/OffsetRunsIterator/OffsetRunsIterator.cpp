#include <hrleFillDomainFromPointList.hpp>
#include <hrleOffsetRunsIterator.hpp>
#include <hrleRunsIterator.hpp>
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

  hrleFillDomainFromPointList(
      data, pointData,
      '.'); // last parameter is the background value to use

  // iterate over hrle structure and output each value and the one above it
  hrleVectorType<hrleIndexType, D> offset(0, -1, 0);
  hrleRunsIterator<hrleDomain<char, D>> it(data);
  hrleOffsetRunsIterator<hrleDomain<char, D>> offsetIt(data, offset);
  while (!it.isFinished()) {
    std::cout << it.getStartIndices() << " " << it.getValue() << "\t\t"
              << offsetIt.getStartIndices() - offset << " "
              << offsetIt.getValue() << std::endl;
    it.next();
    if (compare(it.getStartIndices(), offsetIt.getStartIndices()) > 0)
      offsetIt.next();
  }

  // print data structure in debug output
  // data.print();

  return 0;
}
