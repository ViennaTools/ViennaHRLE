#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <string>

int main() {
  // set dimension for domain
  constexpr int D = 2;

  // set the spacial extension of the domain
  hrleIndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -30;
    max[i] = 30;
  }

  // initialise a domain with extensions min, max in which we can store the type
  // char
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  // the simplest way to fill the domain with data is using a vector of
  // index/value pairs
  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index(0, 0, 0);

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

  // print the data structure to stdout
  data.print(std::cout);

  return 0;
}
