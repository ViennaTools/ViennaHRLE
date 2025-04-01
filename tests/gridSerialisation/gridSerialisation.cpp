#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleTypes.hpp>

#include <vcTestAsserts.hpp>

int main() {
  constexpr int D = 2;
  using namespace viennahrle;

  IndexType min[D] = {-10000000, -100000000};
  IndexType max[D] = {10, 10};
  const double gridDelta = 1.0;
  BoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] =
        BoundaryType::REFLECTIVE_BOUNDARY; // Reflective boundary conditions
  }
  boundaryCons[D - 1] = BoundaryType::INFINITE_BOUNDARY;

  // Initialise example grid
  auto grid = Grid<D>(min, max, gridDelta, boundaryCons);

  // grid.print();
  // std::cout << "\n\n" << std::flush;

  // grid.print();

  // Open file for writing and save serialized level set in it
  std::ofstream fout("test.grid");

  grid.serialize(fout);

  fout.close();

  // Now read the grid in again
  std::ifstream fin("test.grid");

  Grid<D> newGrid;

  newGrid.deserialize(fin);

  fin.close();

  // newGrid.print();

  VC_TEST_ASSERT(grid == newGrid)

  newGrid.print();

  return 0;
}
