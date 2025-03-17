#include <hrleDomainReader.hpp>
#include <hrleDomainWriter.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  using namespace viennahrle;
  omp_set_num_threads(4); // how many threads to use

  // set domain bounds
  IndexType min[D] = {0, -2}, max[D] = {10, 10};

  // declare domain with bounds min,max
  Grid<D> grid(min, max);
  Domain<char, D> data(grid);

  std::vector<std::pair<VectorType<IndexType, D>, char>> pointData;
  VectorType<IndexType, D> index = data.getGrid().getMinIndex();

  std::string dataString = "The quick brown fox jumps over the lazy dog.";

  for (unsigned i = 0; i < dataString.size(); ++i) {
    if (index[0] == data.getGrid().getMaxIndex()[0]) {
      index[0] = data.getGrid().getMinIndex()[0];
      index[1] += 2;
    }
    pointData.emplace_back(index, dataString[i]);
    ++(index[0]);
  }

  FillDomainFromPointList(data, pointData,
                          '#'); // last parameter is the background value to use

  DomainWriter<Domain<char, D>> writer(data);
  writer.setFilePath("test.hrle");
  writer.apply();

  Domain<char, D> newData(grid);
  DomainReader<Domain<char, D>> reader(newData);
  reader.setFilePath("test.hrle");
  reader.apply();

  newData.print();

  return 0;
}
