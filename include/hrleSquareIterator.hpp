#ifndef HRLE_SQUARE_ITERATOR_HPP
#define HRLE_SQUARE_ITERATOR_HPP

#include "hrleOffsetRunsIterator.hpp"
#include "hrleRunsIterator.hpp"

/// This neighbor iterator consists of (2*order+1)^dimension
/// hrleOffsetRunsIterator s for the cartesian neighbors and the center.
/// Whenever one of these iterators reach a defined grid point, the square
/// iterator stops.
template <class hrleDomain> class hrleSquareIterator {

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
  std::vector<hrleOffsetRunsIterator<hrleDomain>> neighborIterators;

  bool isDefined() const {
    if (neighborIterators[centerIndex].isDefined())
      return true;
    for (int i = 0; i < 2 * D; i++) {
      if (neighborIterators[i].isDefined())
        return true;
    }
    return false;
  }

  hrleVectorType<hrleIndexType, D> indexToCoordinate(hrleIndexType index) {
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

  hrleIndexType coordinateToIndex(hrleVectorType<hrleIndexType, D> coordinate) {
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
          hrleOffsetRunsIterator<hrleDomain>(domain, offset, v));
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleSquareIterator<hrleDomain> operator++(int) { // use pre increment instead
    return *this;
  }
  hrleSquareIterator<hrleDomain> operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleSquareIterator(hrleDomain &passedDomain,
                     const hrleVectorType<hrleIndexType, D> &v,
                     const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))) {

    initializeNeigbors(v);
  }

  hrleSquareIterator(hrleDomain &passedDomain, const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleSquareIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleSquareIterator<hrleDomain> &operator--() {
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
  }

  hrleOffsetRunsIterator<hrleDomain> &getNeighbor(int index) {
    return neighborIterators[index];
  }

  hrleOffsetRunsIterator<hrleDomain> &getNeighbor(unsigned index) {
    return neighborIterators[index];
  }

  template <class V>
  hrleOffsetRunsIterator<hrleDomain> &getNeighbor(V &relativeCoordinate) {
    return neighborIterators[coordinateToIndex(relativeCoordinate)];
  }

  hrleOffsetRunsIterator<hrleDomain> &center() {
    return neighborIterators[centerIndex];
  }

  bool isFinished() const { return center().isFinished(); }
};

template <class hrleDomain>
using hrleConstSquareIterator = hrleSquareIterator<const hrleDomain>;

#endif // HRLE_SQUARE_ITERATOR_HPP
