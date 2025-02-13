#include <fstream>
#include <iostream>

#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleGrid.hpp>
#include <hrleTestAsserts.hpp>

constexpr int D = 3;
using DataType = char;

void writeString(std::string fileName, std::string s) {
  std::ofstream out(fileName);
  out << s;
}

void testSerialization(hrleDomain<DataType, D> &domain,
                       bool writeOutput = false) {
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
  hrleGrid<D> newGrid;
  newGrid.deserialize(domainStream);
  // newGrid.deserialize(gridStream);

  hrleDomain<DataType, D> newDomain(&newGrid);
  newDomain.deserialize(domainStream);

  HRLETEST_ASSERT(newGrid == domain.getGrid())
  std::stringstream ss1, ss2;
  domain.print(ss1);
  newDomain.print(ss2);
  HRLETEST_ASSERT(ss1.str() == ss2.str())

  // serialize new grid again to see if all internal variables are also the same
  std::stringstream newDomainStream;
  newGrid.serialize(newDomainStream);
  newDomain.serialize(newDomainStream);

  auto newDomainString = newDomainStream.str();
  HRLETEST_ASSERT(domainString == newDomainString)

  if (writeOutput)
    writeString("newDomainTest.hrle", newDomainString);
}

void fillDomain(hrleDomain<DataType, D> &domain) {
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
      domain, pointData,
      '.'); // last parameter is the background value to use
}

int main() {
  omp_set_num_threads(1);

  std::array<hrleIndexType, D> min = {-20, -20, -20};
  std::array<hrleIndexType, D> max = {20, 20, 20};
  std::array<hrleBoundaryType, D> bounds = {
      hrleBoundaryType::REFLECTIVE_BOUNDARY,
      hrleBoundaryType::INFINITE_BOUNDARY,
      hrleBoundaryType::REFLECTIVE_BOUNDARY};

  hrleGrid<D> grid(min.data(), max.data(), 1.0, bounds.data());

  hrleDomain<DataType, D> domain(&grid);

  // test empty domain
  testSerialization(domain);

  // now fill domain and test again
  fillDomain(domain);
  testSerialization(domain);

  return 0;
}
