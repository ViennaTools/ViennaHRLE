#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseStarIterator.hpp>
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

  // go to the first 'l' and output all neighbors
  constexpr int order = 3;
  hrleConstSparseStarIterator<hrleDomain<char, D>, order> neighborIt(data);
  while (!(neighborIt.getCenter().getValue() == 'l')) {
    neighborIt.next();
  }

  --neighborIt;

  hrleVectorType<hrleIndexType, D> testIndex(-1, 0);
  std::cout << "Value at index " << testIndex << ": "
            << neighborIt.getNeighbor(testIndex).getValue() << std::endl;

  std::cout << "Neighbors at \'l\':" << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << std::string(2 * order, ' ')
              << neighborIt.getNeighbor(2 * D * i + 1).getValue() << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << neighborIt.getNeighbor(2 * D * i + 2).getValue() << " ";
  std::cout << neighborIt.getCenter().getValue() << " ";
  for (unsigned i = 0; i < order; ++i)
    std::cout << neighborIt.getNeighbor(2 * D * i).getValue() << " ";
  std::cout << std::endl;
  for (unsigned i = 0; i < order; ++i)
    std::cout << std::string(2 * order, ' ')
              << neighborIt.getNeighbor(2 * D * i + 3).getValue() << std::endl;

  return 0;
}
