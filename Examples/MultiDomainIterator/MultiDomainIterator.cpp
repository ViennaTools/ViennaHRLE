#include <hrleFillDomainFromPointList.hpp>
#include <hrleSparseMultiIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  omp_set_num_threads(1);

  // set domain bounds
  hrleIndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -5;
    max[i] = 30;
  }

  // declare domain with bounds min,max
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index(0);

  std::string helloString = "Hello, World!";

  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  hrleFillDomainFromPointList(
      data, pointData,
      '.'); // last parameter is the background value to use

  // make second domain to check if multidomainiterator works
  hrleDomain<char, D> data2(grid);

  std::cout << "Filling domain 2" << std::endl;

  pointData.clear();
  index = hrleVectorType<hrleIndexType, D>(3, 1, 0);
  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  hrleFillDomainFromPointList(data2, pointData, '.');

  // iterate over hrle structures and output the values showing the dense data
  // set
  std::cout << "Data Sets filled with background value:" << std::endl;
  hrleConstSparseMultiIterator<hrleDomain<char, D>> it(data);
  it.insertNextDomain(data2);
  int x = data.getGrid().getMinIndex(0);
  int y = data.getGrid().getMinIndex(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      x = data.getGrid().getMinIndex(0);
      std::cout << std::endl;
    }

    for (; x < int(it.getIndex(0)); ++x) {
      std::cout << " ";
    }

    auto definedIts = it.getDefinedIterators();
    for (auto &defined : definedIts) {
      std::cout << defined.second.getValue();
    }

    ++it;
  }

  std::cout << "Finished" << std::endl;

  return 0;
}
