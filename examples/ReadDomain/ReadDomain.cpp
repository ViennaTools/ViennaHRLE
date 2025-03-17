#include <hrleDomain.hpp>
#include <hrleDomainReader.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  using namespace viennahrle;

  Grid<D> grid; // grid will be filled during file read
  Domain<char, D> data(grid);

  DomainReader<Domain<char, D>> reader(data);
  reader.setFilePath("foxDomain.hrle");
  reader.apply();

  data.print();

  // grid now contains the correct values
  std::cout << "Min indices: " << grid.getMinIndex() << std::endl;
  std::cout << "Max Indices: " << grid.getMaxIndex() << std::endl;

  return 0;
}
