#ifndef HRLE_SQUARE_ITERATOR_HPP
#define HRLE_SQUARE_ITERATOR_HPP

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"

/// This neighbor iterator consists of (2*order+1)^dimension
/// hrleSparseOffsetIterator s for the cartesian neighbors and the center.
/// Whenever one of these iterators reach a defined grid point, the square
/// iterator stops.
/// Neighbors are indexed lexicographically from negative cartesian directions:
/// order 1:            order 2:
///                     20 21 22 23 24
/// 6 7 8               15 16 17 18 19
/// 3 4 5               10 11 12 13 14
/// 0 1 2               5  6  7  8  9
///                     0  1  2  3  4
/// center: 4           center: 12
template <class hrleDomain> class hrleSparseBoxIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  const unsigned order;
  const hrleIndexType sideLength;
  const hrleIndexType sliceArea;
  const hrleIndexType centerIndex;
  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleSparseOffsetIterator<hrleDomain>> neighborIterators;

  bool isDefined() const {
    if (neighborIterators[centerIndex].isDefined())
      return true;
    for (int i = 0; i < 2 * D; i++) {
      if (neighborIterators[i].isDefined())
        return true;
    }
    return false;
  }

  hrleVectorType<hrleIndexType, D>
  indexToCoordinate(hrleIndexType index) const {
    hrleVectorType<hrleIndexType, D> coordinate;

    if (D > 2) {
      coordinate[2] = index / sliceArea;
      index = index % sliceArea;
    }

    coordinate[1] = index / sideLength;
    coordinate[0] = index % sideLength;

    // shift to the middle
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] -= order;

    return coordinate;
  }

  template <class V> hrleIndexType coordinateToIndex(V coordinate) const {
    // shift to the middle
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] += order;

    hrleIndexType index = 0;
    if (D > 2) {
      index += coordinate[2] * sliceArea;
    }
    index += coordinate[1] * sideLength;
    index += coordinate[0];

    return index;
  }

  /// push offset iterators lexicographically into std::vector from -order to
  /// +order
  template <class V> void initializeNeigbors(const V &v) {
    const unsigned numNeighbors = unsigned(std::pow((1 + 2 * order), D));

    neighborIterators.reserve(numNeighbors);

    for (unsigned i = 0; i < numNeighbors; ++i) {
      hrleVectorType<hrleIndexType, D> offset = indexToCoordinate(i);
      neighborIterators.push_back(
          hrleSparseOffsetIterator<hrleDomain>(domain, offset, v));
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleSparseBoxIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleSparseBoxIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleSparseBoxIterator(hrleDomain &passedDomain,
                        const hrleVectorType<hrleIndexType, D> &v,
                        const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
        currentCoords(v) {

    initializeNeigbors(v);
  }

  hrleSparseBoxIterator(hrleDomain &passedDomain,
                        const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleSparseBoxIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleSparseBoxIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numNeighbors = int(neighborIterators.size());
    std::vector<bool> increment(numNeighbors + 1, false);
    increment[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> end_coords =
        neighborIterators[centerIndex].getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;

      switch (compare(end_coords, neighborIterators[i].getEndIndices())) {
      case 1:
        end_coords = neighborIterators[i].getEndIndices();
        increment = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        increment[i] = true;
      }
    }

    if (increment[numNeighbors])
      neighborIterators[centerIndex].next();
    for (int i = 0; i < numNeighbors; i++)
      if (increment[i])
        neighborIterators[i].next();

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    const int numNeighbors = neighborIterators.size();
    std::vector<bool> decrement(numNeighbors + 1, false);
    decrement[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> start_coords =
        neighborIterators[centerIndex].getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;
      switch (compare(start_coords, neighborIterators[i].getStartIndices())) {
      case -1:
        start_coords = neighborIterators[i].getStartIndices();
        decrement = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        decrement[i] = true;
      }
    }

    if (decrement[numNeighbors])
      neighborIterators[centerIndex].previous();
    for (int i = 0; i < numNeighbors; i++)
      if (decrement[i])
        neighborIterators[i].previous();

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    return neighborIterators[index];
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    return neighborIterators[index];
  }

  template <class V>
  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(V relativeCoordinate) {
    return neighborIterators[coordinateToIndex(relativeCoordinate)];
  }

  hrleSparseOffsetIterator<hrleDomain> &getCenter() {
    return neighborIterators[centerIndex];
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  unsigned getSize() { return neighborIterators.size(); }

  bool isFinished() const { return getCenter().isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    const unsigned numNeighbors = neighborIterators.size();
    getCenter().goToIndices(v);
    for (int j = 0; j < numNeighbors; ++j) {
      neighborIterators[j].goToIndices(v);
    }
  }

  /// Advances the iterator to position v.
  /// If v is lexicographically higher than the current position
  /// the iterator will be moved back to v.
  /// If v is lexicographically smaller than the current position
  /// then the iterator will be moved until it reaches v
  template <class V> void goToIndicesSequential(const V &v) {
    if (v >= currentCoords) {
      while (v > currentCoords) {
        next();
      }
    } else {
      while (v < currentCoords) {
        previous();
      }
    }
  }
};

template <class hrleDomain>
using hrleConstSparseBoxIterator = hrleSparseBoxIterator<const hrleDomain>;

#endif // HRLE_SQUARE_ITERATOR_HPP
