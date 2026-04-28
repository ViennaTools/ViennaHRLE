#ifndef HRLE_SQUARE_ITERATOR_HPP
#define HRLE_SQUARE_ITERATOR_HPP

#include <array>
#include <utility>

#include "hrleSparseOffsetIterator.hpp"
#include "hrleUtil.hpp"

namespace viennahrle {
using namespace viennacore;
/// This neighbor iterator consists of (2*order+1)^dimension
/// SparseOffsetIterator s for the cartesian neighbors and the center.
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
template <class hrleDomain, int order = 1> class SparseBoxIterator {
public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

private:
  using ValueType = std::conditional_t<std::is_const_v<hrleDomain>,
                                       const typename hrleDomain::ValueType,
                                       typename hrleDomain::ValueType>;
  static constexpr int D = hrleDomain::dimension;
  static constexpr IndexType sideLength = 1 + 2 * order;
  static constexpr IndexType sliceArea = sideLength * sideLength;
  static constexpr auto numNeighbors =
      static_cast<unsigned>(hrleUtil::pow(1 + 2 * order, D));

  hrleDomain &domain;
  const IndexType centerIndex;
  Index<D> currentCoords;
  std::array<OffsetIterator, numNeighbors> neighborIterators;

  static Index<D> indexToCoordinate(IndexType index) {
    Index<D> coordinate;

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

  template <class V> static IndexType coordinateToIndex(V coordinate) {
    // shift to the middle
    for (unsigned i = 0; i < D; ++i)
      coordinate[i] += order;

    IndexType index = 0;
    if (D > 2) {
      index += coordinate[2] * sliceArea;
    }
    index += coordinate[1] * sideLength;
    index += coordinate[0];

    return index;
  }

  template <class V, std::size_t... Is>
  static std::array<OffsetIterator, numNeighbors>
  makeNeighborIteratorsImpl(hrleDomain &passedDomain, const V &v,
                            std::index_sequence<Is...>) {
    return {OffsetIterator(
        passedDomain, indexToCoordinate(static_cast<IndexType>(Is)), v)...};
  }

  template <class V>
  static std::array<OffsetIterator, numNeighbors>
  makeNeighborIterators(hrleDomain &passedDomain, const V &v) {
    return makeNeighborIteratorsImpl(
        passedDomain, v,
        std::make_index_sequence<static_cast<std::size_t>(numNeighbors)>{});
  }

public:
  explicit SparseBoxIterator(hrleDomain &passedDomain)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(domain.getGrid().getMinGridPoint()),
        neighborIterators(makeNeighborIterators(
            passedDomain, passedDomain.getGrid().getMinIndex())) {}

  SparseBoxIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(v),
        neighborIterators(makeNeighborIterators(passedDomain, v)) {}

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  SparseBoxIterator operator++(int) = delete; // use pre increment instead
  SparseBoxIterator operator--(int) = delete; // use pre decrement instead

  SparseBoxIterator &operator++() {
    next();
    return *this;
  }

  SparseBoxIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    std::array<bool, numNeighbors + 1> increment;
    increment.fill(false);
    increment[numNeighbors] = true;

    Index<D> end_coords = neighborIterators[centerIndex].getEndIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;

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
      neighborIterators[centerIndex].next();
    for (int i = 0; i < numNeighbors; i++) {
      if (increment[i])
        neighborIterators[i].next();
    }
    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    std::array<bool, numNeighbors + 1> decrement;
    decrement.fill(false);
    decrement[numNeighbors] = true;

    Index<D> start_coords = neighborIterators[centerIndex].getStartIndices();
    for (int i = 0; i < numNeighbors; i++) {
      if (i == centerIndex)
        continue;
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
      neighborIterators[centerIndex].previous();
    for (int i = 0; i < numNeighbors; i++) {
      if (decrement[i])
        neighborIterators[i].previous();
    }
    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  OffsetIterator &getNeighbor(int index) {
    assert(index >= 0 && index < numNeighbors);
    return neighborIterators[index];
  }

  OffsetIterator const &getNeighbor(int index) const {
    assert(index >= 0 && index < numNeighbors);
    return neighborIterators[index];
  }

  OffsetIterator &getNeighbor(unsigned index) {
    assert(index < numNeighbors);
    return neighborIterators[index];
  }

  OffsetIterator const &getNeighbor(unsigned index) const {
    assert(index < numNeighbors);
    return neighborIterators[index];
  }

  template <class V> OffsetIterator &getNeighbor(V relativeCoordinate) {
    return neighborIterators[coordinateToIndex(relativeCoordinate)];
  }

  OffsetIterator &getCenter() { return neighborIterators[centerIndex]; }

  const OffsetIterator &getCenter() const {
    return neighborIterators[centerIndex];
  }

  const Index<D> &getIndices() { return currentCoords; }

  unsigned getSize() { return neighborIterators.size(); }

  const DomainType &getDomain() { return domain; }

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
      while (v > currentCoords)
        next();
    } else {
      while (v < currentCoords)
        previous();
    }
  }
};

template <class hrleDomain, int order>
using ConstSparseBoxIterator = SparseBoxIterator<const hrleDomain, order>;

} // namespace viennahrle

#endif // HRLE_SQUARE_ITERATOR_HPP
