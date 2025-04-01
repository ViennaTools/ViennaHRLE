#include <fstream>
#include <iostream>

#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>

#include <vcTestAsserts.hpp>

using namespace viennahrle;
constexpr int D = 3;
using DataType = char;

void writeString(std::string const &fileName, std::string const &s) {
  std::ofstream out(fileName);
  out << s;
}

void testSerialization(Domain<DataType, D> &domain, bool writeOutput = false) {
  std::stringstream domainStream;

  // serialise grid
  domain.getGrid().serialize(domainStream);
  // grid.serialize(gridStream);
  domain.serialize(domainStream);

  // std::cout << "Done serializing" << std::endl;
  auto domainString = domainStream.str();
  if (writeOutput)
    writeString("domainTest.hrle", domainString);

  // deserialize
  Grid<D> newGrid;
  newGrid.deserialize(domainStream);
  // newGrid.deserialize(gridStream);

  Domain<DataType, D> newDomain(&newGrid);
  newDomain.deserialize(domainStream);

  VC_TEST_ASSERT(newGrid == domain.getGrid())
  std::stringstream ss1, ss2;
  domain.print(ss1);
  newDomain.print(ss2);
  VC_TEST_ASSERT(ss1.str() == ss2.str())

  // serialize new grid again to see if all internal variables are also the same
  std::stringstream newDomainStream;
  newGrid.serialize(newDomainStream);
  newDomain.serialize(newDomainStream);

  auto newDomainString = newDomainStream.str();
  VC_TEST_ASSERT(domainString == newDomainString)

  if (writeOutput)
    writeString("newDomainTest.hrle", newDomainString);
}

void fillDomain(Domain<DataType, D> &domain) {
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
  for (int i = static_cast<int>(helloString.size()) - 1; i >= 0; --i) {
    pointData.emplace_back(index, helloString[i]);
    index[0] += 1;
    // index[1] += 1;
  }

  FillDomainFromPointList(domain, pointData,
                          '.'); // last parameter is the background value to use
}

int main() {
  omp_set_num_threads(1);

  std::array<IndexType, D> min = {-20, -20, -20};
  std::array<IndexType, D> max = {20, 20, 20};
  std::array<BoundaryType, D> bounds = {BoundaryType::REFLECTIVE_BOUNDARY,
                                        BoundaryType::INFINITE_BOUNDARY,
                                        BoundaryType::REFLECTIVE_BOUNDARY};

  Grid<D> grid(min.data(), max.data(), 1.0, bounds.data());

  Domain<DataType, D> domain(&grid);

  // test empty domain
  testSerialization(domain);

  // now fill domain and test again
  fillDomain(domain);
  testSerialization(domain);

  return 0;
}
