#ifndef HRLE_CROSS_ITERATOR_HPP
#define HRLE_CROSS_ITERATOR_HPP

#include <array>
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
public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

private:
  using ValueType = std::conditional_t<std::is_const_v<hrleDomain>,
                                       const typename hrleDomain::ValueType,
                                       typename hrleDomain::ValueType>;
  static constexpr int D = hrleDomain::dimension;
  static constexpr int numNeighbors = 2 * order * D;

  hrleDomain const &domain;
  Index<D> currentCoords;
  SparseIterator<hrleDomain> centerIterator;
  std::array<OffsetIterator, numNeighbors> neighborIterators;

  static Index<D> makeRelativeIndex(int neighborId) {
    Index<D> relativeIndex(0);
    const int shell = neighborId / (2 * D);
    const int direction = neighborId % (2 * D);

    if (direction < D)
      relativeIndex[direction] = shell + 1;
    else
      relativeIndex[direction - D] = -(shell + 1);

    return relativeIndex;
  }

  template <std::size_t... Is>
  static std::array<OffsetIterator, numNeighbors>
  makeNeighborIteratorsImpl(hrleDomain &passedDomain, const Index<D> &v,
                            std::index_sequence<Is...>) {
    return {OffsetIterator(passedDomain,
                           makeRelativeIndex(static_cast<int>(Is)), v)...};
  }

  static std::array<OffsetIterator, numNeighbors>
  makeNeighborIterators(hrleDomain &passedDomain, const Index<D> &v) {
    return makeNeighborIteratorsImpl(passedDomain, v,
                                     std::make_index_sequence<numNeighbors>{});
  }

public:
  explicit SparseStarIterator(hrleDomain &passedDomain)
      : domain(passedDomain), currentCoords(domain.getGrid().getMinGridPoint()),
        centerIterator(passedDomain),
        neighborIterators(makeNeighborIterators(
            passedDomain, passedDomain.getGrid().getMinIndex())) {
    static_assert(numNeighbors + 1 <= 64,
                  "SparseStarIterator assumes at most 64 iterators");
  }

  SparseStarIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v), centerIterator(passedDomain, v),
        neighborIterators(makeNeighborIterators(passedDomain, v)) {
    static_assert(numNeighbors + 1 <= 64,
                  "SparseStarIterator assumes at most 64 iterators");
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
    Index<D> end_coords = centerIterator.getEndIndices();

    // bit numNeighbors represents the center iterator
    std::uint64_t incrementMask = std::uint64_t{1} << numNeighbors;

    for (int i = 0; i < numNeighbors; ++i) {
      const auto &neighborEnd = neighborIterators[i].getEndIndices();
      const int cmp = Compare(end_coords, neighborEnd);

      if (cmp > 0) {
        end_coords = neighborEnd;
        incrementMask = std::uint64_t{1} << i;
      } else if (cmp == 0) {
        incrementMask |= std::uint64_t{1} << i;
      }
    }

    if (incrementMask & (std::uint64_t{1} << numNeighbors))
      centerIterator.next();

    for (int i = 0; i < numNeighbors; ++i) {
      if (incrementMask & (std::uint64_t{1} << i))
        neighborIterators[i].next();
    }

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    Index<D> start_coords = centerIterator.getStartIndices();

    // bit numNeighbors represents the center iterator
    std::uint64_t decrementMask = std::uint64_t{1} << numNeighbors;

    for (int i = 0; i < numNeighbors; ++i) {
      const auto &neighborStart = neighborIterators[i].getStartIndices();
      const int cmp = Compare(start_coords, neighborStart);

      if (cmp < 0) {
        start_coords = neighborStart;
        decrementMask = std::uint64_t{1} << i;
      } else if (cmp == 0) {
        decrementMask |= std::uint64_t{1} << i;
      }
    }

    if (decrementMask & (std::uint64_t{1} << numNeighbors))
      centerIterator.previous();

    for (int i = 0; i < numNeighbors; ++i) {
      if (decrementMask & (std::uint64_t{1} << i))
        neighborIterators[i].previous();
    }

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  const OffsetIterator &getNeighbor(int index) const {
    return neighborIterators[index];
  }

  OffsetIterator &getNeighbor(int index) {
    return const_cast<OffsetIterator &>(
        const_cast<const SparseStarIterator *>(this)->getNeighbor(index));
  }

  const OffsetIterator &getNeighbor(unsigned index) const {
    return neighborIterators[index];
  }

  OffsetIterator &getNeighbor(unsigned index) {
    return const_cast<OffsetIterator &>(
        const_cast<const SparseStarIterator *>(this)->getNeighbor(index));
  }

  const OffsetIterator &getNeighbor(Index<D> const &relativeIndex) const {
    // check first if it is a valid index
    unsigned char directions = 0;
    unsigned neighborIndex =
        std::numeric_limits<unsigned>::max(); // invalid index
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
    assert(neighborIndex < numNeighbors);

    return neighborIterators[neighborIndex];
  }

  OffsetIterator &getNeighbor(Index<D> const &relativeIndex) {
    return const_cast<OffsetIterator &>(
        const_cast<const SparseStarIterator *>(this)->getNeighbor(
            relativeIndex));
  }

  SparseIterator<hrleDomain> &getCenter() { return centerIterator; }

  const SparseIterator<hrleDomain> &getCenter() const { return centerIterator; }

  const Index<D> &getIndices() const { return currentCoords; }

  const IndexType &getIndices(unsigned i) const { return currentCoords[i]; }

  const DomainType &getDomain() const { return domain; }

  bool isFinished() const { return centerIterator.isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is slower
  /// than goToIndicesSequential for repeated serial calls.
  void goToIndices(const Index<D> &v) {
    centerIterator.goToIndices(v);
    for (auto &it : neighborIterators)
      it.goToIndices(v);
  }

  /// Advances the iterator to position v.
  /// If v is lexicographically higher than the current position
  /// the iterator will be moved back to v.
  /// If v is lexicographically smaller than the current position
  /// then the iterator will be moved until it reaches v
  void goToIndicesSequential(const Index<D> &v) {
    if (v >= currentCoords) {
      while (v > currentCoords)
        next();
    } else {
      while (v < currentCoords)
        previous();
    }
  }
};

template <class hrleDomain, int order>
using ConstSparseStarIterator = SparseStarIterator<const hrleDomain, order>;

} // namespace viennahrle

#endif // HRLE_CROSS_ITERATOR_HPP
