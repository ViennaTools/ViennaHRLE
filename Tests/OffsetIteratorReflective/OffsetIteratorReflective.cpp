#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleTestAsserts.hpp>
#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleDenseIterator.hpp>
#include <hrleSparseBoxIterator.hpp>
#include <hrleSparseIterator.hpp>

constexpr int D = 2;
using DataType = char;
using VectorType = hrleVectorType<hrleIndexType, D>;

const std::array<hrleIndexType, D> min = {-10, -10};
const std::array<hrleIndexType, D> max = {10, 10};

std::array<hrleGrid<D>::boundaryType, D> bounds = {
      hrleGrid<D>::REFLECTIVE_BOUNDARY,
      hrleGrid<D>::INFINITE_BOUNDARY};

void printDenseDomain(hrleDomain<DataType, D> &data) {
  hrleConstDenseIterator<hrleDomain<char, D>> it(data);
  int y = data.getGrid().getMinIndex(1);
  while (!it.isFinished()) {
    if (y < it.getIndex(1)) {
      y = it.getIndex(1);
      std::cout << std::endl;
    }
    std::cout << (it++).getValue() << " ";
  }

  std::cout << std::endl;
}

void fillDomain(hrleDomain<DataType, D> &domain) {
  std::vector<std::pair<hrleVectorType<hrleIndexType, D>, char>> pointData;
  constexpr double radius = 12.;
  constexpr double radius2 = radius * radius;
  constexpr int circleExtent = radius + 2;

  auto &grid = domain.getGrid();
  auto gridMin = grid.getMinGridPoint();
  auto gridMax = grid.getMaxGridPoint();

  hrleVectorType<hrleIndexType, D> index;

  // write an 'x' at all points on the circle
  for(index[1] = -circleExtent; index[1] <= circleExtent; ++index[1]) {
    for(index[0] = gridMin[0]; index[0] <= gridMax[0]; ++index[0]) {
      if(std::abs(double(index[0]*index[0]) + double(index[1]*index[1]) - radius2) < 40.) {
        pointData.push_back(std::make_pair(index, 'x'));
      }
    }
  }

  hrleFillDomainFromPointList(
      domain, pointData,
      '.'); // last parameter is the background value to use
}

bool checkOffset(hrleGrid<D> &grid, VectorType coord, VectorType nbor, VectorType offset) {
  VectorType realOffset;
  auto gridMin = grid.getMinGridPoint();
  auto gridMax = grid.getMaxGridPoint();
  coord += offset;
  for(unsigned i = 0; i < D; ++i) {
    // if neighbour is outside of domain
    if(coord[i] > gridMax[i]) {
      if(bounds[i] == hrleGrid<D>::REFLECTIVE_BOUNDARY) {
        // reflect
        coord[i] = gridMax[i] - (coord[i] - gridMax[i]);
      } else {
        // come in from other side
        coord[i] = gridMin[i] + (coord[i] - gridMax[i] - 1);
      }
    } else if(coord[i] < gridMin[i]) {
      if(bounds[i] == hrleGrid<D>::REFLECTIVE_BOUNDARY) {
        // reflect
        coord[i] = gridMin[i] + (gridMin[i] - coord[i]);
      } else {
        // come in from other side
        coord[i] = gridMax[i] + (coord[i] - gridMin[i] + 1);
      } 
    }

    if(coord[i] != nbor[i]) {
      std::cout << "ERROR: Wrong offset: index=" << coord << "; neighbor=" << nbor << "; correct offset=" << offset[i] << std::endl;
      return false;
    }
  }
  return true;
}

int main() {
  omp_set_num_threads(1);

  hrleGrid<D> grid(min.data(), max.data(), 1.0, bounds.data());

  hrleDomain<DataType, D> domain(&grid);

  fillDomain(domain);

  // march through domain and check if all neighbours are correct
  const unsigned iteratorOrder = 2;
  hrleConstSparseBoxIterator<hrleDomain<char, D>> it(domain, iteratorOrder);
  const unsigned numNeighbors = unsigned(std::pow((1 + 2 * iteratorOrder), D));

  for(; !it.isFinished(); ++it) {
    if(!it.getCenter().isDefined()) {
      continue;
    }
    for(unsigned i = 0; i < numNeighbors; ++i) {
      const auto &neighbor = it.getNeighbor(i);
      if(!neighbor.isDefined()) {
        continue;
      }
      HRLETEST_ASSERT(checkOffset(grid, it.getIndices(), it.getNeighbor(i).getOffsetIndices(), it.getNeighbor(i).getOffset()));
    }
  }

  return 0;
}
