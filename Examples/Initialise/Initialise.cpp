#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <string>

using namespace viennahrle;

int main() {
  // set dimension for domain
  constexpr int D = 2;

  // set the spacial extension of the domain
  IndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -30;
    max[i] = 30;
  }

  // initialise a domain with extensions min, max in which we can store the type
  // char
  Grid<D> grid(min, max);
  Domain<char, D> data(grid);

  // the simplest way to fill the domain with data is using a vector of
  // index/value pairs
  std::vector<std::pair<Index<D>, char>> pointData;
  Index<D> index(0, 0, 0);

  // fill vector with chars from the string "Hello World!" at increasing indices
  std::string helloString = "Hello, World!";
  for (unsigned i = 0; i < helloString.size(); ++i) {
    pointData.push_back(std::make_pair(index, helloString[i]));
    index[0] += 2;
    index[1] += 1;
  }

  // now put the data into an hrleDomain using a helper function
  FillDomainFromPointList(data, pointData,
                          '.'); // last parameter is the background value to use

  // print the data structure to stdout
  data.print(std::cout);

  return 0;
}
