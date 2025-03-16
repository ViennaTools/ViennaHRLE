#pragma once

#include <vcVectorUtil.hpp>

namespace viennahrle {
using namespace viennacore;

// add different implementations for different systems if needed
typedef int IndexType;
typedef unsigned long SizeType;
typedef double CoordType;

template <int D> using Coordinate = VectorType<CoordType, D>;
template <int D> using Index = VectorType<IndexType, D>;
} // namespace viennahrle
