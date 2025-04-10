#ifndef HRLE_SQUARE_ITERATOR_HPP
#define HRLE_SQUARE_ITERATOR_HPP

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

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  static constexpr IndexType sideLength = 1 + 2 * order;
  static constexpr IndexType sliceArea = sideLength * sideLength;
  static constexpr auto numNeighbors =
      static_cast<unsigned>(hrleUtil::pow(1 + 2 * order, D));

  const IndexType centerIndex;
  Index<D> currentCoords;
  std::vector<SparseOffsetIterator<hrleDomain>> neighborIterators;

  Index<D> indexToCoordinate(IndexType index) const {
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

  /// push offset iterators lexicographically into std::vector from -order to
  /// +order
  template <class V> void initializeNeighbors(const V &v) {
    neighborIterators.reserve(numNeighbors);
    for (unsigned i = 0; i < numNeighbors; ++i) {
      auto offset = indexToCoordinate(i);
      neighborIterators.emplace_back(domain, offset, v);
    }
  }

public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

  SparseBoxIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(v) {
    initializeNeighbors(v);
  }

  explicit SparseBoxIterator(hrleDomain &passedDomain)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(domain.getGrid().getMinGridPoint()) {
    initializeNeighbors(passedDomain.getGrid().getMinIndex());
  }

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
    for (int i = 0; i < numNeighbors; i++)
      if (increment[i])
        neighborIterators[i].next();

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
    for (int i = 0; i < numNeighbors; i++)
      if (decrement[i])
        neighborIterators[i].previous();

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  SparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    assert(index >= 0 && index < numNeighbors);
    return neighborIterators[index];
  }

  SparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    assert(index < numNeighbors);
    return neighborIterators[index];
  }

  template <class V>
  SparseOffsetIterator<hrleDomain> &getNeighbor(V relativeCoordinate) {
    return neighborIterators[coordinateToIndex(relativeCoordinate)];
  }

  SparseOffsetIterator<hrleDomain> &getCenter() {
    return neighborIterators[centerIndex];
  }

  const SparseOffsetIterator<hrleDomain> &getCenter() const {
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

template <class hrleDomain, int order>
using ConstSparseBoxIterator = SparseBoxIterator<const hrleDomain, order>;

} // namespace viennahrle

#endif // HRLE_SQUARE_ITERATOR_HPP
