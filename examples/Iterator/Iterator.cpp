#include <hrleCartesianPlaneIterator.hpp>
#include <hrleDenseCellIterator.hpp>
#include <hrleDenseIterator.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseBoxIterator.hpp>
#include <hrleSparseCellIterator.hpp>
#include <hrleSparseStarIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  using namespace viennahrle;

  // set domain bounds
  IndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -10;
    max[i] = 10;
  }

  // declare domain with bounds min,max
  Grid<D> grid(min, max);
  Domain<char, D> data(grid);

  std::vector<std::pair<Index<D>, char>> pointData;
  Index<D> index(0);

  std::string helloString = "Hello, World!";

  for (char &i : helloString) {
    pointData.emplace_back(index, i);
    index[0] += 2;
    index[1] += 1;
  }

  FillDomainFromPointList(data, pointData, '.');
  // last parameter is the background value to use
  // iterate over hrle structure and output the values showing the dense
  // data set
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

  const std::string alphabet = "ABCDEFGHIJK";

  index = Index<D>(0);
  pointData.clear();
  for (int k = 1; k <= alphabet.size(); ++k) {
    pointData.emplace_back(index, alphabet[k - 1]);
    std::cout << index << std::endl;
    index[0] += 1;
    if (k % 4 == 0) {
      index[1] += 1;
      index[0] = 0;
    }
  }

  Domain<char, D> alpha(grid);
  FillDomainFromPointList(alpha, pointData,
                          '.'); // last parameter is the background value to use

  ConstSparseIterator<Domain<char, D>> sparseIt(alpha);
  while (!sparseIt.isFinished()) {
    std::cout << "Idx: " << sparseIt.getStartIndices() << std::endl;
    std::cout << "value: " << sparseIt.getValue() << std::endl;
    sparseIt.next();
  }

  Index<D> idx(0, 0);
  constexpr int order = 2;
  CartesianPlaneIterator<Domain<char, D>, order> it4(alpha, idx);
  std::cout << "Value at " << it4.getCenter().getStartIndices() << ": "
            << it4.getCenter().getValue() << std::endl;
  std::cout << "Neighbors at " << idx << ":" << std::endl;
  const unsigned numNeighbors = it4.getSize();
  std::cout << numNeighbors << std::endl;
  for (unsigned i = 0; i < numNeighbors; ++i)
    std::cout << it4.getNeighbor(i).getValue() << std::endl;

  ConstSparseStarIterator<Domain<char, D>, order> it5(alpha);

  SparseBoxIterator<Domain<char, D>, order> it6(alpha, idx);
  it6.next();

  return 0;
}
