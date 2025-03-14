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
template <class hrleDomain, int order> class hrleSparseStarIterator {

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::hrleValueType,
                             typename hrleDomain::hrleValueType>
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;
  static constexpr int numNeighbors = 2 * order * D;

  hrleDomain &domain;
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

  template <class V> void initializeNeighbors(const V &v) {
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

public:
  using DomainType = hrleDomain;
  using OffsetIterator = hrleSparseOffsetIterator<hrleDomain>;

  hrleSparseStarIterator(hrleDomain &passedDomain,
                         const hrleVectorType<hrleIndexType, D> &v)
      : domain(passedDomain), currentCoords(v),
        centerIterator(passedDomain, v) {

    initializeNeighbors(v);
  }

  explicit hrleSparseStarIterator(hrleDomain &passedDomain)
      : domain(passedDomain), currentCoords(domain.getGrid().getMinGridPoint()),
        centerIterator(passedDomain) {

    initializeNeighbors(passedDomain.getGrid().getMinIndex());
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  hrleSparseStarIterator operator++(int) = delete; // use pre increment instead
  hrleSparseStarIterator operator--(int) = delete; // use pre decrement instead

  hrleSparseStarIterator &operator++() {
    next();
    return *this;
  }

  hrleSparseStarIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    std::array<bool, numNeighbors + 1> increment;
    increment.fill(false);
    increment[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> end_coords =
        centerIterator.getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (
          hrleUtil::Compare(end_coords, neighborIterators[i].getEndIndices())) {
      case 1:
        end_coords = neighborIterators[i].getEndIndices();
        increment.fill(false);
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
    std::array<bool, numNeighbors + 1> decrement;
    decrement.fill(false);
    decrement[numNeighbors] = true;

    hrleVectorType<hrleIndexType, D> start_coords =
        centerIterator.getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (hrleUtil::Compare(start_coords,
                                neighborIterators[i].getStartIndices())) {
      case -1:
        start_coords = neighborIterators[i].getStartIndices();
        decrement.fill(false);
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

  const OffsetIterator &getNeighbor(unsigned index) const {
    return neighborIterators[index];
  }

  OffsetIterator &getNeighbor(unsigned index) {
    return const_cast<OffsetIterator &>(
        const_cast<const hrleSparseStarIterator *>(this)->getNeighbor(index));
  }

  const OffsetIterator &getNeighbor(int index) const {
    return neighborIterators[index];
  }

  OffsetIterator &getNeighbor(const int index) {
    return const_cast<OffsetIterator &>(
        const_cast<const hrleSparseStarIterator *>(this)->getNeighbor(index));
  }

  template <class V> const OffsetIterator &getNeighbor(V &relativeIndex) const {
    // check first if it is a valid index
    unsigned char directions = 0;
    unsigned neighborIndex;
    for (unsigned i = 0; i < D; ++i) {
      assert(abs(relativeIndex[i]) <= order);
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

  template <class V> OffsetIterator &getNeighbor(V &relativeIndex) {
    return const_cast<OffsetIterator &>(
        const_cast<const hrleSparseStarIterator *>(this)->getNeighbor(
            relativeIndex));
  }

  hrleSparseIterator<hrleDomain> &getCenter() { return centerIterator; }

  const hrleSparseIterator<hrleDomain> &getCenter() const {
    return centerIterator;
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() const {
    return currentCoords;
  }

  const DomainType &getDomain() const { return domain; }

  bool isFinished() const { return centerIterator.isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    centerIterator.goToIndices(v);
    for (int i = 0; i < order; ++i) {
      for (int j = 0; j < 2 * D; ++j) {
        neighborIterators[i * 2 * D + j].goToIndices(v);
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

template <class hrleDomain, int order>
using hrleConstSparseStarIterator =
    hrleSparseStarIterator<const hrleDomain, order>;

#endif // HRLE_CROSS_ITERATOR_HPP
