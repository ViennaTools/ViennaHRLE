#include <hrleDomain.hpp>
#include <hrleDomainReader.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;

  // we have to know the type we are reading and the dimension of it
  typedef hrleDomain<char, D> hrleDomainType;

  hrleGrid<D> grid; // grid will be filled during file read
  hrleDomainType data(grid);

  hrleDomainReader<hrleDomainType> reader(data);
  reader.setFilePath("foxDomain.hrle");
  reader.read();

  data.print();

  // grid now contains the correct values
  std::cout << "Min indices: " << grid.getMinIndex() << std::endl;
  std::cout << "Max Indices: " << grid.getMaxIndex() << std::endl;

  return 0;
}
