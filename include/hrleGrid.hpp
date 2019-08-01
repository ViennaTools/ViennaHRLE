#ifndef HRLE_GRID_HPP
#define HRLE_GRID_HPP

#include "hrleCoordType.hpp"
#include "hrleIndexType.hpp"
#include "hrleVectorType.hpp"
#include <cmath>
#include <memory>
#include <vector>

template <int D> class hrleGrid {

public:
  enum boundaryType : unsigned {
    SYMMETRIC_BOUNDARY = 0,
    INFINITE_BOUNDARY = 1,
    PERIODIC_BOUNDARY = 2,
    POS_INFINITE_BOUNDARY = 3,
    NEG_INFINITE_BOUNDARY = 4
  };

  static const int dimension = D;

private:
  static const hrleIndexType INF_EXTENSION =
      std::numeric_limits<hrleIndexType>::max() / 3;

  hrleVectorType<hrleIndexType, D> minIndex, maxIndex,
      indexExtension; // Minimum, maximum, extension of whole grid for all
                      // grid directions

  hrleCoordType gridDelta;

  hrleVectorType<boundaryType, D>
      boundaryConditions; // here the boundary conditions for
                          // all grid directions are stored

  hrleVectorType<hrleIndexType, D> minGridPointCoord,
      maxGridPointCoord; // effective maximum and minimum grid point
                         // coordinates due to periodic boundary conditions the
                         // grid points at opposite boundaries are the same
                         // therefore these indices are less for periodic
                         // boundary conditions for periodic boundary
                         // conditions in direction k,
                         // maxGridPointCoord[k]=max[k]-1 holds

  hrleVectorType<bool, D>
      parities; // here the parities of the grid are stored, the parity is
                // false if the grid_position function is strictly monotonic
                // increasing and true if the grid_position function is
                // strictly monotonic decreasing

public:
  /// empty constructor
  hrleGrid() {
    hrleIndexType min[D] = {}, max[D] = {};
    *this = hrleGrid(min, max);
  }

  hrleGrid(const hrleIndexType *min, const hrleIndexType *max,
           const hrleCoordType delta = 1.0) {
    boundaryType boundaryCons[D];
    for (unsigned i = 0; i < D; ++i) {
      boundaryCons[i] = SYMMETRIC_BOUNDARY; // Reflective boundary conditions
    }

    *this = hrleGrid(min, max, delta, boundaryCons);
  }

  hrleGrid(const hrleIndexType *min, const hrleIndexType *max,
           const hrleCoordType delta, const boundaryType *boundaryCons) {

    gridDelta = delta;

    for (int i = 0; i < D; i++) {
      boundaryConditions[i] = boundaryCons[i];

      minGridPointCoord[i] = min[i];
      maxGridPointCoord[i] = max[i];
      minIndex[i] = min[i];
      maxIndex[i] = max[i];

      if ((boundaryConditions[i] == INFINITE_BOUNDARY) ||
          (boundaryConditions[i] == NEG_INFINITE_BOUNDARY)) {
        minGridPointCoord[i] = -INF_EXTENSION;
        minIndex[i] = -INF_EXTENSION;
      }
      if ((boundaryConditions[i] == INFINITE_BOUNDARY) ||
          (boundaryConditions[i] == POS_INFINITE_BOUNDARY)) {
        maxGridPointCoord[i] = INF_EXTENSION;
        maxIndex[i] = INF_EXTENSION;
      }
      if (boundaryConditions[i] == PERIODIC_BOUNDARY)
        maxGridPointCoord[i]--;
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
  hrleGrid(const hrleGrid &gt)
      : minIndex(gt.minIndex), maxIndex(gt.maxIndex),
        indexExtension(gt.indexExtension),
        boundaryConditions(gt.boundaryConditions),
        minGridPointCoord(gt.minGridPointCoord),
        maxGridPointCoord(gt.maxGridPointCoord), parities(gt.parities) {}

  void print() const {
    std::cout << "minIndex: " << minIndex << std::endl;
    std::cout << "maxIndex: " << maxIndex << std::endl;
    std::cout << "indexExtension: " << indexExtension << std::endl;
    std::cout << "MinGridPoint: " << minGridPointCoord << std::endl;
    std::cout << "MaxGridPoint: " << maxGridPointCoord << std::endl;
    std::cout << "BNC: " << boundaryConditions << std::endl;
  }

  bool parity(int dim) const {
    // parity is false/true if the "grid_position" function in GridTraitsType
    // is monotonic increasing/decreasing respectively for the given grid
    // direction
    return parities[dim];
  }

  bool parity() const {
    // returns the total parity of the grid
    bool b = parity(0);
    for (int i = 1; i < D; i++)
      b ^= parity(i);
    return b;
  }

  hrleIndexType getGridExtent(int dim) const { return indexExtension[dim]; }

  hrleIndexType getMinIndex(int dim) const { return minIndex[dim]; }

  hrleIndexType getMaxIndex(int dim) const { return maxIndex[dim]; }

  inline const hrleVectorType<hrleIndexType, D> &getMinIndex() const {
    return minIndex;
  }

  inline const hrleVectorType<hrleIndexType, D> &getMaxIndex() const {
    return maxIndex;
  }

  /// returns all 2 or 3 boundary conditions
  inline const hrleVectorType<boundaryType, D> &getBoundaryConditions() const {
    return boundaryConditions;
  }

  /// returns the boundary conditions in the specified direction
  inline boundaryType getBoundaryConditions(int dir) const {
    return boundaryConditions[dir];
  }

  /// returns wheter the boundary condition in direction dim is periodic
  bool isBoundaryPeriodic(int dim) const {
    return boundaryConditions[dim] == PERIODIC_BOUNDARY;
  }

  /// returns wheter the boundary condition in direction +dim is infinite
  bool isPosBoundaryInfinite(int dim) const {
    return ((boundaryConditions[dim] == INFINITE_BOUNDARY) ||
            (boundaryConditions[dim] == POS_INFINITE_BOUNDARY));
  }

  /// returns wheter the boundary condition in direction -dim is infinite
  bool isNegBoundaryInfinite(int dim) const {
    return ((boundaryConditions[dim] == INFINITE_BOUNDARY) ||
            (boundaryConditions[dim] == NEG_INFINITE_BOUNDARY));
  }

  /// return whether the point given by vec is within the simulation domain
  template <class V> bool isInDomain(const V &vec) const {
    for (int i = 0; i < D; ++i) {
      if ((vec[i] < getMinIndex(i)) || (vec[i] >= getMaxIndex(i)))
        return false;
    }
    return true;
  }

  /// returns the index[dim] of the maximum grid point which is maxIndex-1 for
  /// periodic BNCs
  inline hrleIndexType getMaxGridPoint(int dim) const {
    return maxGridPointCoord[dim];
  }

  /// returns the index[dim] of the maximum grid point which is minIndex-1 for
  /// periodic BNCs
  inline hrleIndexType getMinGridPoint(int dim) const {
    return minGridPointCoord[dim];
  }

  /// returns the index of the maximum grid point which is maxIndex-1 for
  /// periodic BNCs
  inline const hrleVectorType<hrleIndexType, D> &getMaxGridPoint() const {
    return maxGridPointCoord;
  }

  /// returns the index of the maximum grid point which is minIndex-1 for
  /// periodic BNCs
  inline const hrleVectorType<hrleIndexType, D> &getMinGridPoint() const {
    return minGridPointCoord;
  }

  /// returns whether the point v is at infinity in any dimension
  template <class V> inline bool isAtInfinity(const V &v) const {
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
  inline hrleIndexType localIndex2GlobalIndex(int dim,
                                              hrleIndexType relative_coord,
                                              int cycles,
                                              hrleIndexType offset = 0) const {
    if (cycles == 0) {
      return relative_coord - offset;
    } else {
      if (isBoundaryPeriodic(dim)) {
        return relative_coord - offset +
               cycles * (maxIndex[dim] - minIndex[dim]);
      } else {
        if ((cycles & 1) == 0) { // if cycles is even
          return cycles * indexExtension[dim] + relative_coord - offset;
        } else { // if cycle is odd
          return cycles * indexExtension[dim] + maxIndex[dim] + minIndex[dim] -
                 relative_coord - offset;
        }
      }
    }
  }

  /// This function transforms a global index to the corresponding local
  /// index. The global index is the index of an infinite grid in all
  /// directions. If symmetric or periodic boundary conditions are used, more
  /// than one global index are mapped to the same local index. For the local
  /// index always (getMinGridPoint<=local index<=getMaxGridPoint
  /// holds)
  inline hrleIndexType globalIndex2LocalIndex(int dim,
                                              hrleIndexType absolute_coord,
                                              hrleIndexType offset,
                                              int &cycles) const {
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
    if (((cycles & 1) == 0) || (isBoundaryPeriodic(dim))) {
      return absolute_coord;
    } else {
      return minIndex[dim] + maxIndex[dim] - absolute_coord;
    }
  }

  /// This function transforms a global index to the corresponding local
  /// index. The global index is the index of an infinite grid in all
  /// directions. If symmetric or periodic boundary conditions are used, more
  /// than one global index are mapped to the same local index. For the local
  /// index always (getMinGridPoint<=local index<=getMaxGridPoint
  /// holds)
  inline hrleIndexType globalIndex2LocalIndex(int dim,
                                              hrleIndexType absolute_coord,
                                              hrleIndexType offset = 0) const {
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
    if (b || isBoundaryPeriodic(dim)) {
      return absolute_coord;
    } else {
      return minIndex[dim] + maxIndex[dim] - absolute_coord;
    }
  }

  /// This function transforms a global index vector to the corresponding
  /// local index vector.
  template <class V>
  inline hrleVectorType<hrleIndexType, D>
  globalIndices2LocalIndices(const V &v) const {
    hrleVectorType<hrleIndexType, D> tmp;
    for (int i = 0; i < D; i++) {
      tmp[i] = globalIndex2LocalIndex(i, v[i]);
    }
    return tmp;
  }

  /// Returns the coordinate of the point at index "index" in the direction dir.
  hrleCoordType index2Coordinate(int dir, hrleIndexType Index) const {
    return Index * gridDelta;
  }

  /// This function transforms the coordinate c in respect to the rectilinear
  /// grid into the real coordinates. Non-integer contributions in c are
  /// considered!
  hrleCoordType localIndex2LocalCoordinate(int dir, hrleCoordType c) const {
    hrleCoordType lc = std::floor(c);
    hrleCoordType uc = std::ceil(c);

    if (lc != uc) {
      return (index2Coordinate(dir, static_cast<hrleIndexType>(lc)) *
                  ((uc - c)) + // TODO
              index2Coordinate(dir, static_cast<hrleIndexType>(uc)) *
                  ((c - lc)));
    } else {
      return index2Coordinate(dir, static_cast<hrleIndexType>(lc));
    }
  }

  // This function transforms the coordinate c in respect to the rectilinear
  // grid into the real coordinates.Non-integer contributions in c are
  // considered!
  hrleCoordType globalIndex2GlobalCoordinate(int dir, hrleCoordType c) const {
    hrleCoordType lc = std::floor(c);
    hrleCoordType uc = std::ceil(c);
    if (lc != uc) {
      return (gridPositionOfGlobalIndex(dir, static_cast<hrleIndexType>(lc)) *
                  ((uc - c) / (uc - lc)) + // TODO
              gridPositionOfGlobalIndex(dir, static_cast<hrleIndexType>(uc)) *
                  ((c - lc) / (uc - lc)));
    } else {
      return gridPositionOfGlobalIndex(dir, static_cast<hrleIndexType>(lc));
    }
  }

  /// Transforms a global index vector to a global coordinate vector.
  template <class V>
  hrleVectorType<hrleCoordType, D>
  globalIndices2GlobalCoordinates(const V &v) const {
    hrleVectorType<hrleCoordType, D> tmp;
    for (unsigned i = 0; i < D; ++i)
      tmp[i] = globalIndex2GlobalCoordinate(i, v[i]);
    return tmp;
  }

  /// Transforms a global coordinate in direction dir to a global index.
  hrleIndexType globalCoordinate2GlobalIndex(int dir, hrleCoordType c) const {
    return hrleIndexType(round(c / gridDelta));
  }

  /// Transforms a global coordinate vector to a global index vector.
  template <class V>
  hrleVectorType<hrleIndexType, D>
  globalCoordinates2GlobalIndices(const V &v) const {
    hrleVectorType<hrleIndexType, D> tmp;
    for (unsigned i = 0; i < D; ++i)
      tmp[i] = globalCoordinate2GlobalIndex(i, v[i]);
    return tmp;
  }

  /// Transforms a local coordinate in direction dir to a local index.
  /// Non-integer contributions of c are considered.
  /// a and b allow to restrict the search between two grid indices
  hrleCoordType localCoordinate2LocalIndex(
      int dir, hrleCoordType c, hrleIndexType a,
      hrleIndexType b) const { // TODO global/local coordinates
    hrleCoordType ac = index2Coordinate(dir, a);
    hrleCoordType bc = index2Coordinate(dir, b);

    if (ac > bc) {
      std::swap(ac, bc);
      std::swap(a, b);
    }

    if (c > bc)
      return b;
    if (c < ac)
      return a;

    while (std::abs(a - b) > hrleIndexType(1)) {

      hrleIndexType mid = (a + b) / 2;
      hrleCoordType midc = index2Coordinate(dir, mid);

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

  /// Transforms a global coordinate in direction dir to a global index.
  hrleCoordType localCoordinate2LocalIndex(int dir, hrleCoordType c) const {
    hrleIndexType a = getMinIndex(dir);
    hrleIndexType b = getMaxIndex(dir);
    return localCoordinate2LocalIndex(dir, c, a, b);
  }

  /*hrleCoordType globalCoordinate2GlobalIndex(int dim, hrleCoordType
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
  const hrleCoordType getGridDelta() const { return gridDelta; }

  // Returns the grid position for a global index,
  // taking into account the boundary conditions.
  hrleCoordType
  gridPositionOfGlobalIndex(int dir,
                            hrleIndexType offset) const { // TODO check
    int cycles = 0;

    while (offset < minIndex[dir]) {
      cycles--;
      offset += indexExtension[dir];
    }
    while (offset >= maxIndex[dir]) {
      cycles++;
      offset -= indexExtension[dir];
    }

    if (((cycles & 1) == 0) || isBoundaryPeriodic(dir)) {
      return index2Coordinate(dir, offset) +
             cycles * (index2Coordinate(dir, getMaxIndex(dir)) -
                       index2Coordinate(dir, getMinIndex(dir)));
    } else {
      return (1 + cycles) * index2Coordinate(dir, getMaxIndex(dir)) +
             (1 - cycles) * index2Coordinate(dir, getMinIndex(dir)) -
             index2Coordinate(dir,
                              getMaxIndex(dir) + getMinIndex(dir) - offset);
    }
  }

  // Returns the minimum coordinate of the domain of the grid.
  hrleCoordType getMinLocalCoordinate(int dir) const {
    if (parity(dir)) {
      return index2Coordinate(dir, getMaxIndex(dir));
    } else {
      return index2Coordinate(dir, getMinIndex(dir));
    }
  }

  /// Returns the maximum coordinate of the domain of the grid
  hrleCoordType getMaxLocalCoordinate(int dir) const {
    if (parity(dir)) {
      return index2Coordinate(dir, getMinIndex(dir));
    } else {
      return index2Coordinate(dir, getMaxIndex(dir));
    }
  }

  /// This function increases the index vector v by unity in lexicographical
  /// order.
  template <class V> V incrementIndices(V v) const {
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
  template <class V> V decrementIndices(V v) const {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      if (v[dim] > getMinGridPoint(dim))
        break;
      v[dim] = getMaxGridPoint(dim);
    }
    --v[dim];
    return v;
  }

  // Determine whether index is on border of simulation domain.
  template <class V> bool isBorderPoint(V v) const {
    for (unsigned i = 0; i < D; ++i) {
      if (v[i] <= minIndex[i] || v[i] >= maxIndex[i])
        return true;
    }
    return false;
  }

  // Determine whether point is outside of the domain in direction other than
  // the direction of the infinite boundary.
  template <class V> bool isOutsideOfDomain(V v) const {
    for (unsigned i = 0; i < D; ++i) {
      if (boundaryConditions[i] == INFINITE_BOUNDARY)
        continue;
      if (boundaryConditions[i] == PERIODIC_BOUNDARY && v[i] == maxIndex[i])
        return true;
      if (v[i] < minIndex[i] || v[i] > maxIndex[i])
        return true;
    }
    return false;
  }
};

// };

#endif // HRLE_GRID_HPP
