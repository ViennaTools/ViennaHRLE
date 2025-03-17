#include <hrleCartesianPlaneIterator.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseBoxIterator.hpp>
#include <hrleSparseStarIterator.hpp>
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
    index[0] += 1;
    // index[1] += 1;
  }
  ++index[1];
  index[0] = 0;
  for (int i = int(helloString.size()) - 1; i >= 0; --i) {
    pointData.emplace_back(index, helloString[i]);
    index[0] += 1;
    // index[1] += 1;
  }

  FillDomainFromPointList(data, pointData,
                          '.'); // last parameter is the background value to use

  // go to the first 'l' and output all neighbors
  constexpr int order = 2;
  SparseBoxIterator<Domain<char, D>, order> neighborIt(data);
  while (!(neighborIt.getCenter().getValue() == 'l')) {
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
  ConstSparseStarIterator<Domain<char, D>, order> crossIt(data);
  while (!(crossIt.getCenter().getValue() == 'l')) {
    crossIt.next();
  }

  // Output spatial data for cross iterator
  std::cout << std::endl << std::endl << "CrossIterator: " << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << std::string(2 * order, ' ')
              << crossIt.getNeighbor(2 * D * i + 1).getValue() << std::endl;
  for (int i = order - 1; i >= 0; --i)
    std::cout << crossIt.getNeighbor(2 * D * i + 2).getValue() << " ";
  std::cout << crossIt.getCenter().getValue() << " ";
  for (unsigned i = 0; i < order; ++i)
    std::cout << crossIt.getNeighbor(2 * D * i).getValue() << " ";
  std::cout << std::endl;
  for (unsigned i = 0; i < order; ++i)
    std::cout << std::string(2 * order, ' ')
              << crossIt.getNeighbor(2 * D * i + 3).getValue() << std::endl;

  ConstCartesianPlaneIterator<Domain<char, D>, 1> planeIt(data);

  return 0;
}
