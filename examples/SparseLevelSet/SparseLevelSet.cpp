#include <hrleDenseIterator.hpp>
#include <hrleDomain.hpp>
#include <hrleFillDomainWithSignedDistance.hpp>
#include <iostream>
#include <limits>

// This example shows how to use hrleFillDomainWithSignedDistance
// It was designed to store a signed distance field in the
// hrleDomain, but it can be used to  store any data
// with two different background values

int main() {
  // set dimension for domain
  constexpr int D = 2;
  typedef double valueType;
  using namespace viennahrle;

  // set the spacial extension of the domain
  IndexType min[D], max[D];
  for (unsigned i = 0; i < D; ++i) {
    min[i] = -5;
    max[i] = 5;
  }

  // initialise a domain with extensions min, max in which we can store the type
  // char
  Grid<D> grid(min, max);
  Domain<valueType, D> data(grid);

  // the simplest way to fill the domain with data is using a vector of
  // index/value pairs
  std::vector<std::pair<Index<D>, valueType>> pointData;

  constexpr double radius = 3.3;
  Index<D> centre(0); // initialise with all zeros

  // fill point list with values
  {
    Index<D> index(grid.getMinIndex());
    while (index < grid.getMaxIndex()) {
      // write only first point outside and all inside circle
      VectorType<double, D> distanceVec;
      for (int i = 0; i < D; ++i)
        distanceVec[i] = double(index[i]) - double(centre[i]);
      double distanceFromCentre =
          std::sqrt(DotProduct(distanceVec, distanceVec));
      if (std::abs(distanceFromCentre - radius) < 1.) {
        pointData.emplace_back(index, distanceFromCentre - radius);
      }
      index = grid.incrementIndices(index);
    }
  }

  // now put the data into an hrleDomain using a helper function
  viennahrle::FillDomainWithSignedDistance(data, pointData, -1., 1.);

  // visualise signed distance field on command line
  ConstDenseIterator<Domain<valueType, D>> it(data);
  int y = data.getGrid().getMinIndex(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      std::cout << std::endl;
    }
    std::cout << std::setw(8) << std::setprecision(2) << it.getValue() << " ";
    ++it;
  }

  std::cout << std::endl;

  return 0;
}
