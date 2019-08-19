#ifndef HRLE_CROSS_ITERATOR_HPP
#define HRLE_CROSS_ITERATOR_HPP

#include <cassert>

#include "hrleOffsetRunsIterator.hpp"
#include "hrleRunsIterator.hpp"

/// This neighbor iterator consists of 2*Dimensions hrleOffsetRunsIterator s
/// for the cartesian neighbors and an hrleRunsIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class hrleCrossIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  const unsigned order;
  // hrleVectorType<hrleIndexType, D> startCoords;
  hrleRunsIterator<hrleDomain> centerIterator;
  std::vector<hrleOffsetRunsIterator<hrleDomain>> neighborIterators;

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
            hrleOffsetRunsIterator<hrleDomain>(domain, relativeIndex, v));
      }
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleCrossIterator<hrleDomain> operator++(int) { // use pre increment instead
    return *this;
  }
  hrleCrossIterator<hrleDomain> operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleCrossIterator(hrleDomain &passedDomain,
                    const hrleVectorType<hrleIndexType, D> &v,
                    const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        centerIterator(passedDomain, v) {

    initializeNeigbors(v);
  }

  hrleCrossIterator(hrleDomain &passedDomain, const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder), centerIterator(passedDomain) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleCrossIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleCrossIterator<hrleDomain> &operator--() {
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
  }

  hrleOffsetRunsIterator<hrleDomain> &getNeighbor(int direction) {
    return neighborIterators[direction];
  }

  template <class V>
  hrleOffsetRunsIterator<hrleDomain> &getNeighbor(V &relativeIndex) {
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

  hrleRunsIterator<hrleDomain> &getCenter() { return centerIterator; }

  const hrleVectorType<hrleIndexType, D> &getIndices() {
    return centerIterator.getStartIndices();
  }

  bool isFinished() const { return getCenter().isFinished(); }

  // const hrleVectorType<hrleIndexType, D> &getStartIndices() const {
  //   return startCoords;
  // }
  //
  // hrleIndexType getStartIndex(int dir) const { return startCoords[dir]; }
};

template <class hrleDomain>
using hrleConstCrossIterator = hrleCrossIterator<const hrleDomain>;

#endif // HRLE_CROSS_ITERATOR_HPP
