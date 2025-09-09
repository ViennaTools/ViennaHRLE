#ifndef HRLE_GRID_HPP
#define HRLE_GRID_HPP

#include "hrleTypes.hpp"
#include "hrleUtil.hpp"

#include <cmath>
#include <memory>
#include <ostream>
#include <vector>

namespace viennahrle {
using namespace viennacore;

enum class BoundaryType : unsigned {
  REFLECTIVE_BOUNDARY = 0,
  INFINITE_BOUNDARY = 1,
  PERIODIC_BOUNDARY = 2,
  POS_INFINITE_BOUNDARY = 3,
  NEG_INFINITE_BOUNDARY = 4
};

inline std::ostream &operator<<(std::ostream &os, const BoundaryType &type) {
  switch (type) {
  case BoundaryType::REFLECTIVE_BOUNDARY:
    return os << "REFLECTIVE";
  case BoundaryType::INFINITE_BOUNDARY:
    return os << "INFINITE";
  case BoundaryType::PERIODIC_BOUNDARY:
    return os << "PERIODIC";
  case BoundaryType::POS_INFINITE_BOUNDARY:
    return os << "POS_INFINITE";
  case BoundaryType::NEG_INFINITE_BOUNDARY:
    return os << "NEG_INFINITE";
  default:
    return os << "UNKNOWN";
  }
}

template <int D> class Grid {
public:
  static const int dimension = D;

private:
  static constexpr IndexType INF_EXTENSION =
      std::numeric_limits<IndexType>::max() / 3;

  Index<D> minIndex, maxIndex,
      indexExtension; // Minimum, maximum, extension of whole grid for all
                      // grid directions

  CoordType gridDelta;

  VectorType<BoundaryType, D>
      boundaryConditions; // here the boundary conditions for
                          // all grid directions are stored

  Index<D> minGridPointCoord,
      maxGridPointCoord; // effective maximum and minimum grid point
                         // coordinates due to periodic boundary conditions the
                         // grid points at opposite boundaries are the same
                         // therefore these indices are less for periodic
                         // boundary conditions for periodic boundary
                         // conditions in direction k,
                         // maxGridPointCoord[k]=max[k]-1 holds

  Index<D> minBounds, maxBounds;

  VectorType<bool, D>
      parities; // here the parities of the grid are stored, the parity is
                // false if the grid_position function is strictly monotonic
                // increasing and true if the grid_position function is
                // strictly monotonic decreasing

  static IndexType readGridBoundary(std::istream &stream,
                                    const unsigned gridBoundaryBytes) {
    if (gridBoundaryBytes == 0)
      return INF_EXTENSION;

    char number[sizeof(IndexType)] = {};
    for (unsigned i = 0; i < gridBoundaryBytes; ++i)
      stream.read(number + i, 1);

    if (number[gridBoundaryBytes - 1] >> 7) {
      for (unsigned i = gridBoundaryBytes; i < sizeof(IndexType); ++i)
        number[i] = static_cast<char>(0xFF);
    }
    return *reinterpret_cast<IndexType *>(number);
  }

public:
  /// empty constructor
  Grid() {
    IndexType min[D] = {}, max[D] = {};
    *this = Grid(min, max);
  }

  Grid(const IndexType *min, const IndexType *max,
       const CoordType delta = 1.0) {
    BoundaryType boundaryCons[D];
    for (unsigned i = 0; i < D; ++i)
      boundaryCons[i] = BoundaryType::REFLECTIVE_BOUNDARY;

    *this = Grid(min, max, delta, boundaryCons);
  }

  Grid(const IndexType *min, const IndexType *max, const CoordType delta,
       const BoundaryType *boundaryCons) {

    gridDelta = delta;

    for (int i = 0; i < D; i++) {
      boundaryConditions[i] = boundaryCons[i];

      minGridPointCoord[i] = min[i];
      maxGridPointCoord[i] = max[i];
      minIndex[i] = min[i];
      maxIndex[i] = max[i];
      minBounds[i] = min[i];
      maxBounds[i] = max[i];

      if ((boundaryConditions[i] == BoundaryType::INFINITE_BOUNDARY) ||
          (boundaryConditions[i] == BoundaryType::NEG_INFINITE_BOUNDARY)) {
        minGridPointCoord[i] = -INF_EXTENSION;
        minIndex[i] = -INF_EXTENSION;
      }
      if ((boundaryConditions[i] == BoundaryType::INFINITE_BOUNDARY) ||
          (boundaryConditions[i] == BoundaryType::POS_INFINITE_BOUNDARY)) {
        maxGridPointCoord[i] = INF_EXTENSION;
        maxIndex[i] = INF_EXTENSION;
      }
      if (boundaryConditions[i] == BoundaryType::PERIODIC_BOUNDARY) {
        --maxGridPointCoord[i];
        --maxBounds[i];
      }
    }
    indexExtension = maxIndex - minIndex;

    // initialize parities
    for (int i = 0; i < D; i++) {
      parities[i] =
          gridDelta < 0; // NOTE: maybe in the future different deltas for
                         // each direction or negative deltas are allowed
    }
  }

  /// copy constructor
  // Grid(const Grid &gt)
  //     : minIndex(gt.minIndex), maxIndex(gt.maxIndex),
  //     indexExtension(gt.indexExtension),
  //       gridDelta(gt.gridDelta),
  //       boundaryConditions(gt.boundaryConditions),
  //       minGridPointCoord(gt.minGridPointCoord),
  //       maxGridPointCoord(gt.maxGridPointCoord), parities(gt.parities) {}

  void print() const {
    std::cout << "minIndex: " << minIndex << std::endl;
    std::cout << "maxIndex: " << maxIndex << std::endl;
    std::cout << "indexExtension: " << indexExtension << std::endl;
    std::cout << "MinGridPoint: " << minGridPointCoord << std::endl;
    std::cout << "MaxGridPoint: " << maxGridPointCoord << std::endl;
    std::cout << "BNC: ";
    for (auto &bc : boundaryConditions) {
      std::cout << bc << " ";
    }
    std::cout << std::endl;
    std::cout << "GridDelta: " << gridDelta << std::endl;
  }

  bool parity(int dim) const noexcept {
    // parity is false/true if the "grid_position" function in GridTraitsType
    // is monotonic increasing/decreasing respectively for the given grid
    // direction
    return parities[dim];
  }

  bool parity() const noexcept {
    // returns the total parity of the grid
    bool b = parity(0);
    for (int i = 1; i < D; i++)
      b ^= parity(i);
    return b;
  }

  IndexType getGridExtent(int dim) const noexcept {
    return indexExtension[dim];
  }

  IndexType getMinIndex(int dim) const noexcept { return minIndex[dim]; }

  IndexType getMaxIndex(int dim) const noexcept { return maxIndex[dim]; }

  const Index<D> &getMinIndex() const noexcept { return minIndex; }

  const Index<D> &getMaxIndex() const noexcept { return maxIndex; }

  /// returns all 2 or 3 boundary conditions
  const VectorType<BoundaryType, D> &getBoundaryConditions() const noexcept {
    return boundaryConditions;
  }

  /// returns the boundary conditions in the specified direction
  BoundaryType getBoundaryConditions(int dir) const noexcept {
    return boundaryConditions[dir];
  }

  /// returns whether the boundary condition in direction dim is periodic
  bool isBoundaryPeriodic(int dim) const noexcept {
    return boundaryConditions[dim] == BoundaryType::PERIODIC_BOUNDARY;
  }

  /// returns whether the boundary condition in direction dim is reflective
  bool isBoundaryReflective(int dim) const noexcept {
    return boundaryConditions[dim] == BoundaryType::REFLECTIVE_BOUNDARY;
  }

  /// returns whether the boundary condition in direction +dim is infinite
  bool isPosBoundaryInfinite(int dim) const noexcept {
    return ((boundaryConditions[dim] == BoundaryType::INFINITE_BOUNDARY) ||
            (boundaryConditions[dim] == BoundaryType::POS_INFINITE_BOUNDARY));
  }

  /// returns whether the boundary condition in direction -dim is infinite
  bool isNegBoundaryInfinite(int dim) const noexcept {
    return ((boundaryConditions[dim] == BoundaryType::INFINITE_BOUNDARY) ||
            (boundaryConditions[dim] == BoundaryType::NEG_INFINITE_BOUNDARY));
  }

  /// return whether the point given by vec is within the simulation domain
  template <class V> bool isInDomain(const V &vec) const noexcept {
    for (int i = 0; i < D; ++i) {
      if ((vec[i] < getMinIndex(i)) || (vec[i] >= getMaxIndex(i)))
        return false;
    }
    return true;
  }

  /// returns the index[dim] of the maximum grid point which is maxIndex-1 for
  /// periodic BNCs
  IndexType getMaxGridPoint(int dim) const noexcept {
    return maxGridPointCoord[dim];
  }

  /// returns the index[dim] of the maximum grid point which is minIndex-1 for
  /// periodic BNCs
  IndexType getMinGridPoint(int dim) const noexcept {
    return minGridPointCoord[dim];
  }

  /// returns the index of the maximum grid point which is maxIndex-1 for
  /// periodic BNCs
  const Index<D> &getMaxGridPoint() const noexcept { return maxGridPointCoord; }

  /// returns the index of the maximum grid point which is minIndex-1 for
  /// periodic BNCs
  const Index<D> &getMinGridPoint() const noexcept { return minGridPointCoord; }

  IndexType getMinBounds(int dim) const noexcept { return minBounds[dim]; }

  IndexType getMaxBounds(int dim) const noexcept { return maxBounds[dim]; }

  const Index<D> &getMinBounds() const noexcept { return minBounds; }

  const Index<D> &getMaxBounds() const noexcept { return maxBounds; }

  /// returns whether the point v is at infinity in any dimension
  template <class V> bool isAtInfinity(const V &v) const noexcept {
    for (int i = 0; i < D; i++) {
      if (std::abs(v[i]) == INF_EXTENSION)
        return true;
    }
    return false;
  }

  /// this function transforms a local index to the corresponding global
  /// index. The global index is the index of an infinite grid in all
  /// directions. If symmetric or periodic boundary conditions are used, more
  /// than one global index is mapped to the same local index. For the local
  /// index always (getMinGridPoint<=local index<=getMaxGridPoint
  /// holds)
  IndexType localIndex2GlobalIndex(int dim, IndexType relative_coord,
                                   int cycles, IndexType offset = 0) const {
    if (cycles == 0)
      return relative_coord - offset;

    if (isBoundaryPeriodic(dim))
      return relative_coord - offset + cycles * (maxIndex[dim] - minIndex[dim]);

    if ((cycles & 1) == 0) // if cycles is even
      return cycles * indexExtension[dim] + relative_coord - offset;

    // if cycle is odd
    return cycles * indexExtension[dim] + maxIndex[dim] + minIndex[dim] -
           relative_coord - offset;
  }

  /// This function transforms a global index to the corresponding local
  /// index. The global index is the index of an infinite grid in all
  /// directions. If symmetric or periodic boundary conditions are used, more
  /// than one global index are mapped to the same local index. For the local
  /// index always (getMinGridPoint<=local index<=getMaxGridPoint
  /// holds)
  IndexType globalIndex2LocalIndex(int dim, IndexType absolute_coord,
                                   IndexType offset, int &cycles) const {
    absolute_coord += offset;
    cycles = 0;
    while (absolute_coord < minIndex[dim]) {
      cycles--;
      absolute_coord += indexExtension[dim];
    }
    while (absolute_coord >= maxIndex[dim]) {
      cycles++;
      absolute_coord -= indexExtension[dim];
    }

    if ((cycles & 1) == 0 || isBoundaryPeriodic(dim))
      return absolute_coord;

    return minIndex[dim] + maxIndex[dim] - absolute_coord;
  }

  /// This function transforms a global index to the corresponding local
  /// index. The global index is the index of an infinite grid in all
  /// directions. If symmetric or periodic boundary conditions are used, more
  /// than one global index are mapped to the same local index. For the local
  /// index always (getMinGridPoint<=local index<=getMaxGridPoint
  /// holds)
  IndexType globalIndex2LocalIndex(int dim, IndexType absolute_coord,
                                   IndexType offset = 0) const {
    absolute_coord += offset;
    bool b = true;
    while (absolute_coord < minIndex[dim]) {
      b = !b;
      absolute_coord += indexExtension[dim];
    }
    while (absolute_coord >= maxIndex[dim]) {
      b = !b;
      absolute_coord -= indexExtension[dim];
    }

    if (b || isBoundaryPeriodic(dim))
      return absolute_coord;

    return minIndex[dim] + maxIndex[dim] - absolute_coord;
  }

  /// This function transforms a global index vector to the corresponding
  /// local index vector.
  template <class V> Index<D> globalIndices2LocalIndices(const V &v) const {
    Index<D> tmp;
    for (int i = 0; i < D; i++)
      tmp[i] = globalIndex2LocalIndex(i, v[i]);
    return tmp;
  }

  /// Returns the coordinate of the point at index "index" in the direction
  /// dir.
  CoordType index2Coordinate(IndexType idx) const noexcept {
    return idx * gridDelta;
  }

  /// This function transforms the coordinate c in respect to the rectilinear
  /// grid into the real coordinates. Non-integer contributions in c are
  /// considered!
  CoordType localIndex2LocalCoordinate(CoordType c) const {
    const CoordType lc = std::floor(c);
    const CoordType uc = std::ceil(c);

    if (lc != uc)
      return index2Coordinate(static_cast<IndexType>(lc)) * (uc - c) +
             index2Coordinate(static_cast<IndexType>(uc)) * (c - lc);

    return index2Coordinate(static_cast<IndexType>(lc));
  }

  /// This function transforms the coordinate c in respect to the rectilinear
  /// grid into the real coordinates.Non-integer contributions in c are
  /// considered!
  CoordType globalIndex2GlobalCoordinate(int dir, CoordType c) const {
    CoordType lc = std::floor(c);
    CoordType uc = std::ceil(c);
    if (lc != uc) {
      return gridPositionOfGlobalIndex(dir, static_cast<IndexType>(lc)) *
                 ((uc - c) / (uc - lc)) +
             gridPositionOfGlobalIndex(dir, static_cast<IndexType>(uc)) *
                 ((c - lc) / (uc - lc));
    } else {
      return gridPositionOfGlobalIndex(dir, static_cast<IndexType>(lc));
    }
  }

  /// Transforms a global index vector to a global coordinate vector.
  template <class V>
  VectorType<CoordType, D> globalIndices2GlobalCoordinates(const V &v) const {
    VectorType<CoordType, D> tmp;
    for (unsigned i = 0; i < D; ++i)
      tmp[i] = globalIndex2GlobalCoordinate(i, v[i]);
    return tmp;
  }

  /// Transforms a global coordinate in direction dir to a global index.
  IndexType globalCoordinate2GlobalIndex(const CoordType c) const noexcept {
    return static_cast<IndexType>(round(c / gridDelta));
  }

  /// Transforms a global coordinate vector to a global index vector.
  template <class V>
  Index<D> globalCoordinates2GlobalIndices(const V &v) const {
    Index<D> tmp;
    for (unsigned i = 0; i < D; ++i)
      tmp[i] = globalCoordinate2GlobalIndex(v[i]);
    return tmp;
  }

  /// Transforms a local coordinate in direction dir to a local index.
  /// Non-integer contributions of c are considered.
  /// a and b allow to restrict the search between two grid indices
  CoordType localCoordinate2LocalIndex(
      CoordType c, IndexType a,
      IndexType b) const { // TODO global/local coordinates
    CoordType ac = index2Coordinate(a);
    CoordType bc = index2Coordinate(b);

    if (ac > bc) {
      std::swap(ac, bc);
      std::swap(a, b);
    }

    if (c > bc)
      return b;
    if (c < ac)
      return a;

    while (std::abs(a - b) > 1) {

      const IndexType mid = (a + b) / 2;
      const CoordType midc = index2Coordinate(mid);

      if (c <= midc) {
        b = mid;
        bc = midc;
      } else {
        a = mid;
        ac = midc;
      }
    }

    return a + (b - a) * ((c - ac) / (bc - ac));
  }

  /// Transforms a local coordinate in direction dir to a local index.
  CoordType localCoordinate2LocalIndex(int dir, CoordType c) const {
    IndexType a = getMinIndex(dir);
    IndexType b = getMaxIndex(dir);
    return localCoordinate2LocalIndex(c, a, b);
  }

  /// Transforms a global coordinate in direction dir to a global index.
  CoordType globalCoordinate2LocalIndex(int dir, CoordType c) const {
    // a is closest grid point less than c
    IndexType a = std::floor(c / gridDelta);
    // get normalised distance within grid
    const CoordType l = c - (a * gridDelta);
    a = globalIndex2LocalIndex(dir, a);

    return static_cast<CoordType>(a) + l / gridDelta;
  }

  /*CoordType globalCoordinate2GlobalIndex(int dim, CoordType
  absolute_coord) const { int cycles=0; if (isBoundaryPeriodic(dim)) { while
  (absolute_coord<getMinLocalCoordinate(dim)) { cycles--;
              absolute_coord+=(getMaxLocalCoordinate(dim)-getMinLocalCoordinate(dim));
          }
          while (absolute_coord>=getMaxLocalCoordinate(dim)) {
              cycles++;
              absolute_coord-=(getMaxLocalCoordinate(dim)-getMinLocalCoordinate(dim));
          }
          if (parity(dim)) cycles=-cycles;
          return localCoordinate2LocalIndex(dim,
  absolute_coord)+cycles*indexExtension[dim]; } else { do { if
  (absolute_coord<getMinLocalCoordinate(dim)) { cycles--;
                  absolute_coord=getMinLocalCoordinate(dim)+(getMinLocalCoordinate(dim)-absolute_coord);
                  continue;
              }
              if (absolute_coord>getMaxLocalCoordinate(dim)) {
                  cycles++;
                  absolute_coord=getMaxLocalCoordinate(dim)+(getMaxLocalCoordinate(dim)-absolute_coord);
                  continue;
              }
          } while(false);
      }
      if (parity(dim)) cycles=-cycles;

      if ((cycles & 1)==0) {
          return localCoordinate2LocalIndex(dim,
  absolute_coord)+cycles*indexExtension[dim]; } else { return
  cycles*indexExtension[dim]+maxIndex[dim]+minIndex[dim]-localCoordinate2LocalIndex(dim,
  absolute_coord);
      }
  }*/

  // const GridTraitsType &grid_traits() const { return GridTraits; }

  /// returns the grid point separation
  CoordType getGridDelta() const noexcept { return gridDelta; }

  /// Returns the grid position for a global index,
  /// taking into account the boundary conditions.
  CoordType gridPositionOfGlobalIndex(int dir,
                                      IndexType offset) const { // TODO check
    int cycles = 0;

    while (offset < minIndex[dir]) {
      cycles--;
      offset += indexExtension[dir];
    }
    while (offset >= maxIndex[dir]) {
      cycles++;
      offset -= indexExtension[dir];
    }

    if ((cycles & 1) == 0 || isBoundaryPeriodic(dir)) {
      return index2Coordinate(offset) +
             cycles * (index2Coordinate(getMaxIndex(dir)) -
                       index2Coordinate(getMinIndex(dir)));
    } else {
      return (1 + cycles) * index2Coordinate(getMaxIndex(dir)) +
             (1 - cycles) * index2Coordinate(getMinIndex(dir)) -
             index2Coordinate(getMaxIndex(dir) + getMinIndex(dir) - offset);
    }
  }

  // Returns the minimum coordinate of the domain of the grid.
  CoordType getMinLocalCoordinate(int dir) const {
    if (parity(dir)) {
      return index2Coordinate(getMaxIndex(dir));
    } else {
      return index2Coordinate(getMinIndex(dir));
    }
  }

  /// Returns the maximum coordinate of the domain of the grid
  CoordType getMaxLocalCoordinate(int dir) const {
    if (parity(dir)) {
      return index2Coordinate(getMinIndex(dir));
    } else {
      return index2Coordinate(getMaxIndex(dir));
    }
  }

  /// This function increases the index vector v by unity in lexicographical
  /// order.
  template <class V> V incrementIndices(V v) const noexcept {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (v[dim] < getMaxGridPoint(dim))
        break;
      v[dim] = getMinGridPoint(dim);
    }
    ++v[dim];
    return v;
  }

  /// This function decreases the index vector v by unity in lexicographical
  /// order.
  template <class V> V decrementIndices(V v) const noexcept {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (v[dim] > getMinGridPoint(dim))
        break;
      v[dim] = getMaxGridPoint(dim);
    }
    --v[dim];
    return v;
  }

  /// Determine whether index is on border of simulation domain.
  template <class V> bool isBorderPoint(V v) const noexcept {
    for (unsigned i = 0; i < D; ++i) {
      if (v[i] <= minIndex[i] || v[i] >= maxIndex[i])
        return true;
    }
    return false;
  }

  /// Determine whether point is outside the domain in direction other than
  /// the direction of the infinite boundary.
  template <class V> bool isOutsideOfDomain(V v) const noexcept {
    for (unsigned i = 0; i < D; ++i) {
      if (boundaryConditions[i] == BoundaryType::INFINITE_BOUNDARY)
        continue;
      if (boundaryConditions[i] == BoundaryType::PERIODIC_BOUNDARY &&
          v[i] == maxIndex[i])
        return true;
      if (v[i] < minIndex[i] || v[i] > maxIndex[i])
        return true;
    }
    return false;
  }

  /// Serialize the grid
  std::ostream &serialize(std::ostream &stream) const {
    using namespace hrleUtil;
    // set identifier
    stream << "hrleGrid";
    // GRID PROPERTIES
    IndexType bounds[2 * D];
    // get the bounds to save
    for (unsigned i = 0; i < D; ++i) {
      if (getBoundaryConditions(i) == BoundaryType::INFINITE_BOUNDARY ||
          getBoundaryConditions(i) == BoundaryType::NEG_INFINITE_BOUNDARY) {
        // set to zero as it will be changed during deserialization
        bounds[2 * i] = 0;
      } else {
        bounds[2 * i] = minIndex[i];
      }

      if (getBoundaryConditions(i) == BoundaryType::INFINITE_BOUNDARY ||
          getBoundaryConditions(i) == BoundaryType::POS_INFINITE_BOUNDARY) {
        bounds[2 * i + 1] = 0;
      } else {
        bounds[2 * i + 1] = maxIndex[i];
      }
    }

    {
      // find number of bytes needed to represent the highest grid extent
      char gridBoundaryBytes = 0;
      for (unsigned i = 0; i < 2 * D; ++i) {
        if (bounds[i] != 0) {
          gridBoundaryBytes =
              std::max(getByteSizeOfNumber(bounds[i]), gridBoundaryBytes);
        }
      }
      gridBoundaryBytes =
          std::min(gridBoundaryBytes, static_cast<char>(sizeof(IndexType)));

      // grid properties
      stream.write(reinterpret_cast<char *>(&gridBoundaryBytes), sizeof(char));
      for (int dim = D - 1; dim >= 0; --dim) {
        stream.write((char *)&bounds[2 * dim], gridBoundaryBytes);
        stream.write((char *)&bounds[2 * dim + 1], gridBoundaryBytes);
        char boundaryCondition = (char)getBoundaryConditions(dim);
        stream.write((char *)&boundaryCondition, 1);
      }
      stream.write((char *)&gridDelta, sizeof(double));
    }

    return stream;
  }

  /// Deserialize from serialized input
  std::istream &deserialize(std::istream &stream) {
    // check identifier
    char identifier[9] = {}; // 1 more for string constructor
    stream.read(identifier, 8);
    if (std::string(identifier).compare(0, 8, "hrleGrid")) {
      std::cout << "Reading Grid from stream failed. Header could not be found."
                << std::endl;
      return stream;
    }

    char gridBoundaryBytes = 0;
    stream.read(&gridBoundaryBytes, 1);
    // READ GRID
    IndexType gridMin[D], gridMax[D];
    BoundaryType boundaryCons[D];
    double gridDelta = 0.;
    for (int i = D - 1; i >= 0; --i) {
      gridMin[i] = readGridBoundary(stream, gridBoundaryBytes);
      gridMax[i] = readGridBoundary(stream, gridBoundaryBytes);
      char condition = 0;
      stream.read(&condition, 1);
      boundaryCons[i] = static_cast<BoundaryType>(condition);
    }
    stream.read((char *)&gridDelta, sizeof(double));
    // initialize new grid
    *this = Grid(gridMin, gridMax, gridDelta, boundaryCons);

    return stream;
  }

  bool operator==(const Grid &other) const {
    for (unsigned i = 0; i < D; ++i) {
      if (minIndex[i] != other.minIndex[i] ||
          maxIndex[i] != other.maxIndex[i] ||
          boundaryConditions[i] != other.boundaryConditions[i] ||
          minGridPointCoord[i] != other.minGridPointCoord[i] ||
          maxGridPointCoord[i] != other.maxGridPointCoord[i])
        return false;
    }
    return true;
  }

  bool operator!=(const Grid &other) const { return !(*this == other); }
};
} // namespace viennahrle

#endif // HRLE_GRID_HPP
