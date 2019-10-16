#ifndef HRLE_CROSS_ITERATOR_HPP
#define HRLE_CROSS_ITERATOR_HPP

#include <cassert>

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"

/// This neighbor iterator consists of 2*Dimensions hrleSparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class hrleSparseStarIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  const unsigned order;
  hrleVectorType<hrleIndexType, D> currentCoords;
  hrleSparseIterator<hrleDomain> centerIterator;
  std::vector<hrleSparseOffsetIterator<hrleDomain>> neighborIterators;

  bool isDefined() const {
    if (centerIterator.isDefined())
      return true;
    for (int i = 0; i < 2 * D; i++) {
      if (neighborIterators[i].isDefined())
        return true;
    }
    return false;
  }

  template <class V> void initializeNeigbors(const V &v) {
    for (int i = 0; i < int(order); ++i) {
      for (int j = 0; j < 2 * D; ++j) {
        hrleVectorType<hrleIndexType, D> relativeIndex(hrleIndexType(0));
        if (j < D)
          relativeIndex[j] = i + 1;
        else
          relativeIndex[j - D] = -(i + 1);
        neighborIterators.push_back(
            hrleSparseOffsetIterator<hrleDomain>(domain, relativeIndex, v));
      }
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleSparseStarIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleSparseStarIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleSparseStarIterator(hrleDomain &passedDomain,
                         const hrleVectorType<hrleIndexType, D> &v,
                         const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder), currentCoords(v),
        centerIterator(passedDomain, v) {

    initializeNeigbors(v);
  }

  hrleSparseStarIterator(hrleDomain &passedDomain,
                         const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        currentCoords(domain.getGrid().getMinGridPoint()),
        centerIterator(passedDomain) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleSparseStarIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleSparseStarIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numNeighbors = 2 * order * D;
    std::vector<bool> increment(numNeighbors + 1, false);
    increment[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> end_coords =
        centerIterator.getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (compare(end_coords, neighborIterators[i].getEndIndices())) {
      case 1:
        end_coords = neighborIterators[i].getEndIndices();
        increment = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        increment[i] = true;
      }
    }

    if (increment[numNeighbors])
      centerIterator.next();
    for (int i = 0; i < numNeighbors; ++i)
      if (increment[i])
        neighborIterators[i].next();

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    const int numNeighbors = 2 * order * D;
    std::vector<bool> decrement(numNeighbors + 1, false);
    decrement[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> start_coords =
        centerIterator.getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (compare(start_coords, neighborIterators[i].getStartIndices())) {
      case -1:
        start_coords = neighborIterators[i].getStartIndices();
        decrement = std::vector<bool>(numNeighbors + 1, false);
      case 0:
        decrement[i] = true;
      }
    }

    if (decrement[numNeighbors])
      centerIterator.previous();
    for (int i = 0; i < numNeighbors; ++i) {
      if (decrement[i])
        neighborIterators[i].previous();
    }
    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    return neighborIterators[index];
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    return neighborIterators[index];
  }

  template <class V>
  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(V &relativeIndex) {
    // check first if it is a valid index
    unsigned char directions = 0;
    unsigned neighborIndex;
    for (unsigned i = 0; i < D; ++i) {
      assert(unsigned(abs(relativeIndex[i])) <= order);
      if (relativeIndex[i] != 0) {
        ++directions;
        if (relativeIndex[i] > 0)
          neighborIndex = 2 * D * (relativeIndex[i] - 1) + i;
        else
          neighborIndex = 2 * D * ((-relativeIndex[i]) - 1) + D + i;
      }
    }
    assert(directions == 1);

    return neighborIterators[neighborIndex];
  }

  hrleSparseIterator<hrleDomain> &getCenter() { return centerIterator; }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  bool isFinished() const { return centerIterator.isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is be slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    centerIterator.goToIndices(v);
    for (int i = 0; i < int(order); ++i) {
      for (int j = 0; j < 2 * D; ++j) {
        neighborIterators[j].goToIndices(v);
      }
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

  // const hrleVectorType<hrleIndexType, D> &getStartIndices() const {
  //   return startCoords;
  // }
  //
  // hrleIndexType getStartIndex(int dir) const { return startCoords[dir]; }
};

template <class hrleDomain>
using hrleConstSparseStarIterator = hrleSparseStarIterator<const hrleDomain>;

#endif // HRLE_CROSS_ITERATOR_HPP
