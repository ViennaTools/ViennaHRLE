#include <hrleDomainReader.hpp>
#include <hrleDomainWriter.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
#ifdef _OPENMP
  omp_set_num_threads(4); // how many threads to use
#endif

  typedef hrleDomain<char, D> hrleDomainType;

  // set domain bounds
  hrleIndexType min[D] = {0, -2}, max[D] = {10, 10};

  // declare domain with bounds min,max
  hrleGrid<D> grid(min, max);
  hrleDomainType data(grid);

  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index = data.getGrid().getMinIndex();

  std::string dataString = "The quick brown fox jumps over the lazy dog.";

  for (unsigned i = 0; i < dataString.size(); ++i) {
    if (index[0] == data.getGrid().getMaxIndex()[0]) {
      index[0] = data.getGrid().getMinIndex()[0];
      index[1] += 2;
    }
    pointData.push_back(std::make_pair(index, dataString[i]));
    ++(index[0]);
  }

  hrleFillDomainFromPointList(
      data, pointData,
      '#'); // last parameter is the background value to use

  hrleDomainWriter<hrleDomainType> writer(data);
  writer.setFilePath("test.hrle");
  writer.apply();

  hrleDomainType newData(grid);
  hrleDomainReader<hrleDomainType> reader(newData);
  reader.setFilePath("test.hrle");
  reader.apply();

  newData.print();

  return 0;
}
