#include <hrleFillDomainFromPointList.hpp>
#include <hrleIterator.hpp>
#include <hrleRunsIterator.hpp>
#include <iostream>
#include <string>

int main() {

  constexpr int D = 2;
  omp_set_num_threads(4); // how many threads to use

  // set domain bounds
  hrleIndexType min[D] = {0, -2}, max[D] = {10, 10};

  // declare domain with bounds min,max
  hrleGrid<D> grid(min, max);
  hrleDomain<char, D> data(grid);

  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  hrleVectorType<hrleIndexType, D> index = data.getGrid().getMinIndex();

  std::string dataString = "The quick brown fox jumps over the lazy dog.";

  for (hrleIndexType i = 0; i < dataString.size(); ++i) {
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

  // split the datastructure across domainSegments for easier parallelization
  data.segment();

// replace all 'o' and 'e' characters with the ID of the thread which changed it
#pragma omp parallel
  {
    int p = 0;
#ifdef _OPENMP
    p = omp_get_thread_num();
#endif
    // define iterator starting at the corresponding domainSegment
    hrleRunsIterator<hrleDomain<char, D>> it(
        data, (p == 0) ? grid.getMinIndex() : data.getSegmentation()[p - 1]);
    // store the end of the domainSegment
    hrleVectorType<hrleIndexType, D> endOfSegment =
        (p != static_cast<int>(data.getSegmentation().size()))
            ? data.getSegmentation()[p]
            : grid.getMaxIndex();

    for (; it.getStartIndices() < endOfSegment; it.next()) {
      if (it.getValue() == 'o' || it.getValue() == 'e')
        it.getValue() = 48 + p; // character for thread_num
    }
  }

  hrleConstIterator<hrleDomain<char, D>> pit(data);
  int y = data.getGrid().getMinIndex(1);
  while (!pit.isFinished()) {
    if (y < pit.getIndex(1)) {
      y = pit.getIndex(1);
      std::cout << std::endl;
    }
    // post in/decrement is available for simple iterators
    std::cout << (pit++).getValue() << " ";
  }

  std::cout << pit.getValue() << std::endl;

  return 0;
}
