#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleTestAsserts.hpp>

int main() {
  constexpr int D = 2;

  hrleIndexType min[D] = {-10000000, -100000000};
  hrleIndexType max[D] = {10, 10};
  const double gridDelta = 1.0;
  hrleBoundaryType boundaryCons[D];
  for (unsigned i = 0; i < D; ++i) {
    boundaryCons[i] =
        hrleBoundaryType::REFLECTIVE_BOUNDARY; // Reflective boundary conditions
  }
  boundaryCons[D - 1] = hrleBoundaryType::INFINITE_BOUNDARY;

  // Initialise example grid
  auto grid = hrleGrid<D>(min, max, gridDelta, boundaryCons);

  // grid.print();
  // std::cout << "\n\n" << std::flush;

  // grid.print();

  // Open file for writing and save serialized level set in it
  std::ofstream fout("test.grid");

  grid.serialize(fout);

  fout.close();

  // Now read the grid in again
  std::ifstream fin("test.grid");

  hrleGrid<D> newGrid;

  newGrid.deserialize(fin);

  fin.close();

  // newGrid.print();

  HRLETEST_ASSERT(grid == newGrid)

  return 0;
}
