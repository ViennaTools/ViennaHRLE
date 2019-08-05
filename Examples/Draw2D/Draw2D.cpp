#include <hrleFillDomainFromPointList.hpp>
#include <hrleRunsIterator.hpp>

#include <iostream>
#include <string>

// only output each run once, to draw the sparse data
template <class hrleDomain> void draw2D(hrleDomain &domain) {
  // iterate over hrle structure and output the corresponding value
  hrleRunsIterator<hrleDomain> it(domain, true);
  int y = domain.getGrid().getMaxIndex(1);
  int x = domain.getGrid().getMinIndex(0);
  std::cout << it.getValue();
  while (!it.isFinished()) {
    if (y > it.getStartIndices()[1]) {
      x = domain.getGrid().getMinIndex(0);
      for (; y > it.getStartIndices()[1]; --y)
        std::cout << std::endl;
    }

    for (; x < it.getStartIndices()[0]; ++x)
      std::cout << " ";

    std::cout << it.getValue();
    --it;
  }
  std::cout << std::endl << std::endl;
}

int main() {
  // set dimension for domain
  constexpr int D = 2;

  // set the spacial extension of the domain
  hrleIndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -15;
    max[i] = 15;
  }

  // initialise a domain with extensions min, max in which we can store the type
  // char
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  // the simplest way to fill the domain with data is using a vector of
  // index/value pairs
  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index(-10, -10, 0);

  // fill vector with chars from the string "Hello World!" at increasing indices
  std::string helloString = "Hello, World!";
  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  // now put the data into an hrleDomain using a helper function
  hrleFillDomainFromPointList(
      data, pointData,
      '.'); // last parameter is the background value to use

  // plot the 2D data in the terminal
  draw2D(data);

  return 0;
}
