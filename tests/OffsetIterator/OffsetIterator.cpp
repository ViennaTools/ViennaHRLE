#include <fstream>
#include <iostream>

#include <hrleDenseIterator.hpp>
#include <hrleDomain.hpp>
#include <hrleFillDomainFromPointList.hpp>
#include <hrleGrid.hpp>
#include <hrleSparseBoxIterator.hpp>
#include <hrleSparseIterator.hpp>

#include <vcTestAsserts.hpp>

using namespace viennahrle;
constexpr int D = 2;
using DataType = char;

void printDenseDomain(Domain<DataType, D> &data) {
  ConstDenseIterator<Domain<char, D>> it(data);
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

void fillDomain(Domain<DataType, D> &domain) {
  std::vector<std::pair<Index<D>, char>> pointData;
  constexpr double radius = 12.;
  constexpr double radius2 = radius * radius;
  constexpr int circleExtent = radius + 2;

  auto &grid = domain.getGrid();
  auto gridMin = grid.getMinGridPoint();
  auto gridMax = grid.getMaxGridPoint();

  Index<D> index;

  // write an 'x' at all points on the circle
  for (index[1] = -circleExtent; index[1] <= circleExtent; ++index[1]) {
    for (index[0] = gridMin[0]; index[0] <= gridMax[0]; ++index[0]) {
      if (std::abs(double(index[0] * index[0]) + double(index[1] * index[1]) -
                   radius2) < 40.) {
        pointData.push_back(std::make_pair(index, 'x'));
      }
    }
  }

  FillDomainFromPointList(domain, pointData,
                          '.'); // last parameter is the background value to use
}

bool checkOffset(const Grid<D> &grid, Index<D> coord, Index<D> nbor,
                 Index<D> offset) {
  Index<D> realOffset;
  auto gridMin = grid.getMinGridPoint();
  auto gridMax = grid.getMaxGridPoint();
  coord += offset;
  for (int i = 0; i < D; ++i) {
    // if neighbour is outside of domain
    if (coord[i] > gridMax[i]) {
      if (grid.isBoundaryReflective(i)) {
        // reflect
        coord[i] = gridMax[i] - (coord[i] - gridMax[i]);
      } else {
        // come in from other side
        coord[i] = gridMin[i] + (coord[i] - gridMax[i] - 1);
      }
    } else if (coord[i] < gridMin[i]) {
      if (grid.isBoundaryReflective(i)) {
        // reflect
        coord[i] = gridMin[i] + (gridMin[i] - coord[i]);
      } else {
        // come in from other side
        coord[i] = gridMax[i] + (coord[i] - gridMin[i] + 1);
      }
    }

    if (coord[i] != nbor[i]) {
      std::cout << "ERROR: Wrong offset: index=" << coord
                << "; neighbor=" << nbor << "; correct offset=" << offset[i]
                << std::endl;
      return false;
    }
  }
  return true;
}

void runTest(Domain<DataType, D> &domain) {
  auto &grid = domain.getGrid();

  fillDomain(domain);

  // march through domain and check if all neighbors are correct
  constexpr int iteratorOrder = 2;
  ConstSparseBoxIterator<Domain<char, D>, iteratorOrder> it(domain);
  constexpr int numNeighbors = hrleUtil::pow((1 + 2 * iteratorOrder), D);

  for (; !it.isFinished(); ++it) {
    if (!it.getCenter().isDefined()) {
      continue;
    }
    for (int i = 0; i < numNeighbors; ++i) {
      const auto &neighbor = it.getNeighbor(i);
      if (!neighbor.isDefined()) {
        continue;
      }
      VC_TEST_ASSERT(checkOffset(grid, it.getIndices(),
                                 it.getNeighbor(i).getOffsetIndices(),
                                 it.getNeighbor(i).getOffset()));
    }
  }
}

int main() {
  omp_set_num_threads(1);

  const std::array<IndexType, D> min = {-10, -10};
  const std::array<IndexType, D> max = {10, 10};

  // Test for periodic boundary condition
  {
    std::array<BoundaryType, D> bounds = {BoundaryType::PERIODIC_BOUNDARY,
                                          BoundaryType::INFINITE_BOUNDARY};

    Grid<D> grid(min.data(), max.data(), 1.0, bounds.data());
    Domain<DataType, D> domain(&grid);

    runTest(domain);
    std::cout << "Passed periodic" << std::endl;
  }

  // Test for reflective boundary condition
  {
    std::array<BoundaryType, D> bounds = {BoundaryType::REFLECTIVE_BOUNDARY,
                                          BoundaryType::INFINITE_BOUNDARY};

    Grid<D> grid(min.data(), max.data(), 1.0, bounds.data());
    Domain<DataType, D> domain(&grid);

    runTest(domain);
    std::cout << "Passed reflective" << std::endl;
  }

  return 0;
}
