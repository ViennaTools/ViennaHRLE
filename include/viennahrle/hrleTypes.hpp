#pragma once

#include <vcVectorUtil.hpp>

namespace viennahrle {
using namespace viennacore;

// add different implementations for different systems if needed
typedef int IndexType;
typedef unsigned long SizeType;
typedef double CoordType;

using Coordinate = VectorType<CoordType, 3>;
template <size_t D> using Index = VectorType<IndexType, D>;
} // namespace viennahrle
