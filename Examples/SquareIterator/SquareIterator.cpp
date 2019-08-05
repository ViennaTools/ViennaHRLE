#include <hrleCrossIterator.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleSquareIterator.hpp>
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
  const int order = 2;
  hrleSquareIterator<hrleDomain<char, D>> neighborIt(data, order);
  while (!(neighborIt.center().getValue() == 'l')) {
    neighborIt.next();
  }

  // output spatial data for square iterator
  std::cout << "Neighbors at \'l\':" << std::endl;
  std::cout << "Square Iterator: " << std::endl;
  for (int j = (1 + 2 * order) * 2 * order; j >= 0; j -= (1 + 2 * order)) {
    for (int i = 0; i < (1 + 2 * order); ++i) {
      // std::cout << i+j;
      std::cout << neighborIt.getNeighbor(i + j).getValue() << " "
                << std::flush;
    }
    std::cout << std::endl;
  }

  // make cross iterator for comparison and advance to the first 'l'
  hrleConstCrossIterator<hrleDomain<char, D>> crossIt(data, order);
  while (!(crossIt.center().getValue() == 'l')) {
    crossIt.next();
  }

  // Output spatial data for cross iterator
  std::cout << std::endl << std::endl << "CrossIterator: " << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << std::string(2 * order, ' ')
              << crossIt.getNeighbor(2 * D * i + 1).getValue() << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << crossIt.getNeighbor(2 * D * i + 2).getValue() << " ";
  std::cout << crossIt.center().getValue() << " ";
  for (unsigned i = 0; i < order; ++i)
    std::cout << crossIt.getNeighbor(2 * D * i).getValue() << " ";
  std::cout << std::endl;
  for (unsigned i = 0; i < order; ++i)
    std::cout << std::string(2 * order, ' ')
              << crossIt.getNeighbor(2 * D * i + 3).getValue() << std::endl;

  return 0;
}
