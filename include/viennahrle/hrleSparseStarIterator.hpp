#ifndef HRLE_CROSS_ITERATOR_HPP
#define HRLE_CROSS_ITERATOR_HPP

#include <cassert>

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"
#include "hrleTypes.hpp"

namespace viennahrle {
using namespace viennacore;
/// This neighbor iterator consists of 2*Dimensions SparseOffsetIterator s
/// for the cartesian neighbors and an SparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain, int order> class SparseStarIterator {

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  static constexpr int D = hrleDomain::dimension;
  static constexpr int numNeighbors = 2 * order * D;

  hrleDomain &domain;
  Index<D> currentCoords;
  SparseIterator<hrleDomain> centerIterator;
  std::vector<SparseOffsetIterator<hrleDomain>> neighborIterators;

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
    neighborIterators.reserve(numNeighbors);
    for (int i = 0; i < order; ++i) {
      for (int j = 0; j < 2 * D; ++j) {
        Index<D> relativeIndex(0);
        if (j < D)
          relativeIndex[j] = i + 1;
        else
          relativeIndex[j - D] = -(i + 1);
        neighborIterators.emplace_back(domain, relativeIndex, v);
      }
    }
  }

public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

  SparseStarIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v),
        centerIterator(passedDomain, v) {
    initializeNeighbors(v);
  }

  explicit SparseStarIterator(hrleDomain &passedDomain)
      : domain(passedDomain), currentCoords(domain.getGrid().getMinGridPoint()),
        centerIterator(passedDomain) {
    initializeNeighbors(passedDomain.getGrid().getMinIndex());
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  SparseStarIterator operator++(int) = delete; // use pre increment instead
  SparseStarIterator operator--(int) = delete; // use pre decrement instead

  SparseStarIterator &operator++() {
    next();
    return *this;
  }

  SparseStarIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    std::array<bool, numNeighbors + 1> increment;
    increment.fill(false);
    increment[numNeighbors] = true;

    Index<D> end_coords = centerIterator.getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (Compare(end_coords, neighborIterators[i].getEndIndices())) {
      case 1:
        end_coords = neighborIterators[i].getEndIndices();
        increment.fill(false);
      case 0:
        increment[i] = true;
      default:
        break;
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

    Index<D> start_coords = centerIterator.getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      switch (Compare(start_coords, neighborIterators[i].getStartIndices())) {
      case -1:
        start_coords = neighborIterators[i].getStartIndices();
        decrement.fill(false);
      case 0:
        decrement[i] = true;
      default:
        break;
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
        const_cast<const SparseStarIterator *>(this)->getNeighbor(index));
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
        const_cast<const SparseStarIterator *>(this)->getNeighbor(
            relativeIndex));
  }

  SparseIterator<hrleDomain> &getCenter() { return centerIterator; }

  const SparseIterator<hrleDomain> &getCenter() const { return centerIterator; }

  const Index<D> &getIndices() const { return currentCoords; }

  const DomainType &getDomain() const { return domain; }

  bool isFinished() const { return centerIterator.isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    centerIterator.goToIndices(v);
    for (auto &neighbor : neighborIterators)
      neighbor.goToIndices(v);
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

  // const Index<D> &getStartIndices() const {
  //   return startCoords;
  // }
  //
  // IndexType getStartIndex(int dir) const { return startCoords[dir]; }
};

template <class hrleDomain, int order>
using ConstSparseStarIterator = SparseStarIterator<const hrleDomain, order>;

} // namespace viennahrle

#endif // HRLE_CROSS_ITERATOR_HPP
