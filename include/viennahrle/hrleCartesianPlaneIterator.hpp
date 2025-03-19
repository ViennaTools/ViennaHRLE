#ifndef HRLE_CARTESIAN_PLANE_ITERATOR_HPP
#define HRLE_CARTESIAN_PLANE_ITERATOR_HPP

#include "hrleSparseOffsetIterator.hpp"
#include "hrleUtil.hpp"

#include <memory>

/// This neighbor iterator is an adaptation of the box iterator, which only
/// moves the iterators that lie in a cartesian plane (2-dimensional linear
/// subspace created by all pairs of basis vectors). Whenever one of these
/// iterators reach a defined grid point, the iterator stops. In 2D this
/// iterator is equivalent to a hrleBoxIterator.
///
/// The indices of the planeCoords coords array for D > 2 are like star
/// iterators along the x-axis except for the slice that contains the center
/// iterator which contains the whole plane
///
///   5        12 13 14        19
/// 2 3 4       9 10 11     16 17 18
///   1         6  7  8        15
///            center: 10

namespace viennahrle {
using namespace viennacore;
template <class hrleDomain, int order = 1> class CartesianPlaneIterator {

  static constexpr int D = hrleDomain::dimension;
  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  hrleDomain &domain;
  static constexpr IndexType sideLength = 1 + 2 * order;
  static constexpr IndexType sliceArea = sideLength * sideLength;
  static constexpr auto numNeighbors =
      static_cast<unsigned>(hrleUtil::pow(1 + 2 * order, D));
  static constexpr auto numPlaneCoords = static_cast<unsigned>(
      hrleUtil::pow(1 + 2 * order, D) - 8 * hrleUtil::pow(order, D));

  const IndexType centerIndex;
  Index<D> currentCoords;
  std::vector<std::unique_ptr<SparseOffsetIterator<hrleDomain>>>
      neighborIterators;
  std::vector<unsigned> planeCoords;

  bool isDefined() const {
    if (neighborIterators[centerIndex].isDefined())
      return true;
    for (int i = 0; i < 2 * D; i++) {
      if (neighborIterators[i].isDefined())
        return true;
    }
    return false;
  }

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
    // TODO: consider empty indices
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

  // CURRENT PLANE ITERATOR IMPLEMENTATION:
  // Create a Box iterator
  // Create a planeCords array that contains the indices of the cartesian planes
  // Only increment the iterators with indices in the planeCoords array
  template <class V> void initializeNeighbors(const V &v) {
    neighborIterators.reserve(numNeighbors);
    planeCoords.reserve(numPlaneCoords);
    for (unsigned i = 0; i < numNeighbors; ++i) {
      Index<D> offset = indexToCoordinate(i);
      // get coordinates that only lie in planes
      int coords = 0;
      for (int j = 0; j < D; j++) {
        if (offset[j] != 0)
          coords++;
      }
      // if more than 2 coordinates are not equal to 0 the point lies not in a
      // cartesian plane
      if (coords < 3) {
        neighborIterators.push_back(
            std::make_unique<SparseOffsetIterator<hrleDomain>>(domain, offset,
                                                               v));
        planeCoords.push_back(i);
      } else {
        neighborIterators.push_back(nullptr);
      }
    }
  }

public:
  using DomainType = hrleDomain;

  CartesianPlaneIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(v) {
    initializeNeighbors(v);
  }

  explicit CartesianPlaneIterator(hrleDomain &passedDomain)
      : domain(passedDomain), centerIndex(coordinateToIndex(Index<D>(0))),
        currentCoords(domain.getGrid().getMinGridPoint()) {
    initializeNeighbors(passedDomain.getGrid().getMinIndex());
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  CartesianPlaneIterator operator++(int) = delete; // use pre increment instead
  CartesianPlaneIterator operator--(int) = delete; // use pre decrement instead

  CartesianPlaneIterator &operator++() {
    next();
    return *this;
  }

  CartesianPlaneIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    const auto numPlaneNeighbours = planeCoords.size();
    std::vector<bool> increment(numPlaneNeighbours + 1, false);
    increment[numPlaneNeighbours] = true;

    Index<D> end_coords = neighborIterators[centerIndex]->getEndIndices();

    // iterate over plane coords
    for (size_t i = 0; i < numPlaneNeighbours; i++) {
      if (planeCoords[i] == centerIndex)
        continue;

      switch (Compare(end_coords,
                      neighborIterators[planeCoords[i]]->getEndIndices())) {
      case 1:
        end_coords = neighborIterators[planeCoords[i]]->getEndIndices();
        increment = std::vector<bool>(numPlaneNeighbours + 1, false);
      case 0:
        increment[i] = true;
      default:
        break;
      }
    }

    if (increment[numPlaneNeighbours])
      neighborIterators[centerIndex]->next();

    for (size_t i = 0; i < numPlaneNeighbours; i++)
      if (increment[i])
        neighborIterators[planeCoords[i]]->next();

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    const auto numPlaneNeighbours = planeCoords.size();
    std::vector<bool> decrement(numPlaneNeighbours + 1, false);
    decrement[numPlaneNeighbours] = true;

    Index<D> start_coords = neighborIterators[centerIndex]->getStartIndices();
    for (size_t i = 0; i < numPlaneNeighbours; i++) {
      if (i == centerIndex)
        continue;

      switch (Compare(start_coords,
                      neighborIterators[planeCoords[i]]->getStartIndices())) {
      case -1:
        start_coords = neighborIterators[planeCoords[i]]->getStartIndices();
        decrement = std::vector<bool>(numPlaneNeighbours + 1, false);
      case 0:
        decrement[i] = true;
      default:
        break;
      }
    }

    if (decrement[numPlaneNeighbours])
      neighborIterators[centerIndex]->previous();
    for (size_t i = 0; i < numPlaneNeighbours; i++)
      if (decrement[i])
        neighborIterators[planeCoords[i]]->previous();

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  SparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    return *(neighborIterators[planeCoords[index]]);
  }

  SparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    return *(neighborIterators[planeCoords[index]]);
  }

  template <class V>
  SparseOffsetIterator<hrleDomain> &getNeighbor(V relativeCoordinate) {
    return *(neighborIterators[coordinateToIndex(relativeCoordinate)]);
  }

  SparseOffsetIterator<hrleDomain> &getCenter() {
    return *(neighborIterators[centerIndex]);
  }

  const Index<D> &getIndices() { return currentCoords; }

  unsigned getSize() const { return planeCoords.size(); }

  const DomainType &getDomain() { return domain; }

  bool isFinished() const { return getCenter().isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    const unsigned numPlaneNeighbours = planeCoords.size();
    getCenter().goToIndices(v);
    for (int j = 0; j < numPlaneNeighbours; ++j) {
      neighborIterators[planeCoords[j]]->goToIndices(v);
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
using ConstCartesianPlaneIterator =
    CartesianPlaneIterator<const hrleDomain, order>;

} // namespace viennahrle

#endif // HRLE_CARTESIAN_PLANE_ITERATOR_HPP