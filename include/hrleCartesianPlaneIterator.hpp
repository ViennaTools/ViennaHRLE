#ifndef HRLE_CARTESIAN_PLANE_ITERATOR_HPP
#define HRLE_CARTESIAN_PLANE_ITERATOR_HPP

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"

/// This neighbor iterator is an adabtation of the box iterator, which only
/// moves the iterators that lie in a cartesian plan (2 dimensional linear
/// subspace created by all pairs of basis vectors). Whenever one of these
/// iterators reach a defined grid point, the iterator stops. In 2D this
/// iterator is equivalent to a hrleBoxIterator.
///
/// The indices of the planeCoords coords array for D > 2 are like star
/// iterators along the x axis exept for the slice that contains the center
/// iterator which contains the whole plane
///
///   5        12 13 14        19
/// 2 3 4       9 10 11     16 17 18
///   1         6  7  8        15
///            center: 10

template <class hrleDomain> class hrleCartesianPlaneIterator {

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
  std::vector<std::unique_ptr<hrleSparseOffsetIterator<hrleDomain>>>
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
    // TODO: Beachte empty indices
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

  // CURRENT PLANE ITERATOR IMPLEMENTATION:
  // Create a Box iterator
  // Create a planeCords array that contains the indices of the cartesian planes
  // Only increment the iterators with indices in the planeCoords array
  template <class V> void initializeNeigbors(const V &v) {
    const unsigned numNeighbors = unsigned(std::pow((1 + 2 * order), D));

    const unsigned numPlaneCoords =
        unsigned(std::pow((1 + 2 * order), D) - (8 * std::pow(order, D)));

    neighborIterators.reserve(numNeighbors);

    planeCoords.reserve(numPlaneCoords);

    for (unsigned i = 0; i < numNeighbors; ++i) {
      hrleVectorType<hrleIndexType, D> offset = indexToCoordinate(i);
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
            std::unique_ptr<hrleSparseOffsetIterator<hrleDomain>>(
                new hrleSparseOffsetIterator<hrleDomain>(domain, offset, v)));
        planeCoords.push_back(i);
      } else {
        neighborIterators.push_back(nullptr);
      }
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleCartesianPlaneIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleCartesianPlaneIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleCartesianPlaneIterator(hrleDomain &passedDomain,
                             const hrleVectorType<hrleIndexType, D> &v,
                             const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
        currentCoords(v) {

    initializeNeigbors(v);
  }

  hrleCartesianPlaneIterator(hrleDomain &passedDomain,
                             const unsigned passedOrder = 1)
      : domain(passedDomain), order(passedOrder),
        sideLength(1 + 2 * passedOrder), sliceArea(sideLength * sideLength),
        centerIndex(coordinateToIndex(hrleVectorType<hrleIndexType, D>(0))),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    initializeNeigbors(passedDomain.getGrid().getMinIndex());
  }

  hrleCartesianPlaneIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleCartesianPlaneIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numPlaneNeighbours = int(planeCoords.size());
    std::vector<bool> increment(numPlaneNeighbours + 1, false);
    increment[numPlaneNeighbours] = true;

    hrleVectorType<hrleIndexType, D> end_coords =
        neighborIterators[centerIndex]->getEndIndices();

    // itearte over plane coords
    for (int i = 0; i < numPlaneNeighbours; i++) {
      if (planeCoords[i] == centerIndex)
        continue;

      switch (compare(end_coords,
                      neighborIterators[planeCoords[i]]->getEndIndices())) {
      case 1:
        end_coords = neighborIterators[planeCoords[i]]->getEndIndices();
        increment = std::vector<bool>(numPlaneNeighbours + 1, false);
      case 0:
        increment[i] = true;
      }
    }

    if (increment[numPlaneNeighbours])
      neighborIterators[centerIndex]->next();

    for (int i = 0; i < numPlaneNeighbours; i++)
      if (increment[i])
        neighborIterators[planeCoords[i]]->next();

    currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    const int numPlaneNeighbours = int(planeCoords.size());
    std::vector<bool> decrement(numPlaneNeighbours + 1, false);
    decrement[numPlaneNeighbours] = true;

    hrleVectorType<hrleIndexType, D> start_coords =
        neighborIterators[centerIndex]->getStartIndices();
    for (int i = 0; i < numPlaneNeighbours; i++) {
      if (i == centerIndex)
        continue;

      switch (compare(start_coords,
                      neighborIterators[planeCoords[i]]->getStartIndices())) {
      case -1:
        start_coords = neighborIterators[planeCoords[i]]->getStartIndices();
        decrement = std::vector<bool>(numPlaneNeighbours + 1, false);
      case 0:
        decrement[i] = true;
      }
    }

    if (decrement[numPlaneNeighbours])
      neighborIterators[centerIndex]->previous();
    for (int i = 0; i < numPlaneNeighbours; i++)
      if (decrement[i])
        neighborIterators[planeCoords[i]]->previous();

    currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(int index) {
    return *(neighborIterators[planeCoords[index]]);
  }

  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(unsigned index) {
    return *(neighborIterators[planeCoords[index]]);
  }

  template <class V>
  hrleSparseOffsetIterator<hrleDomain> &getNeighbor(V relativeCoordinate) {
    return *(neighborIterators[coordinateToIndex(relativeCoordinate)]);
  }

  hrleSparseOffsetIterator<hrleDomain> &getCenter() {
    return *(neighborIterators[centerIndex]);
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  unsigned getSize() { return planeCoords.size(); }

  bool isFinished() const { return getCenter().isFinished(); }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is be slower
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

template <class hrleDomain>
using hrleConstCartesianPlaneIterator =
    hrleCartesianPlaneIterator<const hrleDomain>;

#endif // HRLE_CARTESIAN_PLANE_ITERATOR_HPP