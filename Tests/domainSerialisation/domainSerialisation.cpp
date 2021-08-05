#include <fstream>
#include <iostream>

#include <hrleGrid.hpp>
#include <hrleDomain.hpp>
#include <hrleTestAsserts.hpp>

int main() {
    std::array<hrleIndexType, 3> min = {-5, -5, -5};
    std::array<hrleIndexType, 3> max = {5, 5, 5};
    std::array<hrleGrid<3>::boundaryType, 3> bounds = {hrleGrid<3>::REFLECTIVE_BOUNDARY,
                                                       hrleGrid<3>::INFINITE_BOUNDARY,
                                                       hrleGrid<3>::REFLECTIVE_BOUNDARY};
    std::stringstream ss;

    {
        hrleGrid<3> grid(min.data(), max.data(), 1.0, bounds.data());
        hrleDomain<double, 3> domain(&grid);
        domain.serialize(ss);
    }

    // deserialize
    {
        hrleGrid<3> grid(min.data(), max.data(), 4.0, bounds.data());
        hrleDomain<double, 3> domain(&grid);
        domain.deserialize(ss);
        domain.print();
    }

    return 0;
}
