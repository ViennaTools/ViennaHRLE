#ifndef HRLE_DENSE_CELL_ITERATOR_HPP
#define HRLE_DENSE_CELL_ITERATOR_HPP

#include <cassert>

#include "hrleSparseOffsetIterator.hpp"

/// This neighbor iterator consists of 2*Dimensions hrleSparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class hrleDenseCellIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleSparseOffsetIterator<hrleDomain>> cornerIterators;
  hrleVectorType<hrleIndexType, D> minIndex, maxIndex;

  template <class V> void initialize(const V &v) {
    for (int i = 0; i < (1 << D); ++i) {
      cornerIterators.push_back(hrleSparseOffsetIterator<hrleDomain>(
          domain, BitMaskToVector<D, hrleIndexType>(i), v));
    }
  }

  void incrementIndices(hrleVectorType<hrleIndexType, D> &v) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      bool posInfinite = domain.getGrid().isPosBoundaryInfinite(dim);
      bool negInfinite = domain.getGrid().isNegBoundaryInfinite(dim);
      if (v[dim] < (posInfinite ? domain.getMaxRunBreak(dim)
                                : domain.getGrid().getMaxGridPoint(dim)))
        break;
      v[dim] = (negInfinite ? domain.getMinRunBreak(dim) - 3
                            : domain.getGrid().getMinGridPoint(dim));
    }
    ++v[dim];
  }

  void decrementIndices(hrleVectorType<hrleIndexType, D> &v) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      bool posInfinite = domain.getGrid().isPosBoundaryInfinite(dim);
      bool negInfinite = domain.getGrid().isNegBoundaryInfinite(dim);
      if (v[dim] > (negInfinite ? domain.getMinRunBreak(dim)
                                : domain.getGrid().getMinGridPoint(dim)))
        break;
      v[dim] = (posInfinite ? domain.getMaxRunBreak(dim) + 3
                            : domain.getGrid().getMaxGridPoint(dim));
    }
    --v[dim];
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleDenseCellIterator<hrleDomain>
  operator++(int); // use pre increment instead
  hrleDenseCellIterator<hrleDomain>
  operator--(int); // use pre decrement instead

public:
  hrleDenseCellIterator(hrleDomain &passedDomain,
                        const hrleVectorType<hrleIndexType, D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(currentCoords);
    if (!isDefined())
      next();
  }

  hrleDenseCellIterator(hrleDomain &passedDomain, bool reverse = false)
      : domain(passedDomain),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    auto &grid = domain.getGrid();
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = (grid.isNegBoundaryInfinite(i)) ? domain.getMinRunBreak(i)
                                                    : grid.getMinBounds(i);
      maxIndex[i] = (grid.isPosBoundaryInfinite(i)) ? domain.getMaxRunBreak(i)
                                                    : grid.getMaxBounds(i);
    }

    if (reverse)
      currentCoords = maxIndex;
    else
      currentCoords = minIndex;

    initialize(currentCoords);
    if (!isDefined())
      next();
  }

  bool isDefined() const {
    for (unsigned i = 0; i < D; ++i) {
      if (!domain.getGrid().isBoundaryPeriodic(i) &&
          currentCoords[i] == domain.getGrid().getMaxGridPoint(i)) {
        return false;
      }
    }
    for (int i = 0; i < 2 * D; i++) {
      if (cornerIterators[i].isDefined())
        return true;
    }
    return false;
  }

  hrleDenseCellIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleDenseCellIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numCorners = 1 << D;
    std::vector<bool> increment(numCorners, false);
    increment[0] = true;

    hrleVectorType<hrleIndexType, D> end_coords = currentCoords;
    // cornerIterators[0].getEndIndices();
    for (int i = 0; i < numCorners; i++) {
      switch (compare(end_coords, cornerIterators[i].getEndIndices())) {
      case 1:
        // end_coords = cornerIterators[i].getEndIndices();
        // increment = std::vector<bool>(numCorners, false);
      case 0:
        cornerIterators[i].next();
      }
    }

    incrementIndices(currentCoords);
  }

  void previous() {
    const int numCorners = 1 << D;
    std::vector<bool> decrement(numCorners, false);
    decrement[0] = true;

    hrleVectorType<hrleIndexType, D> start_coords = currentCoords;
    // cornerIterators[0].getStartIndices();
    for (int i = 0; i < numCorners; i++) {
      switch (compare(start_coords, cornerIterators[i].getStartIndices())) {
      case -1:
        // start_coords = cornerIterators[i].getStartIndices();
        // decrement = std::vector<bool>(numCorners, false);
      case 0:
        cornerIterators[i].previous();
      }
    }

    decrementIndices(currentCoords);
  }

  hrleSparseOffsetIterator<hrleDomain> &getCorner(unsigned index) {
    return cornerIterators[index];
  }

  hrleSparseOffsetIterator<hrleDomain> &getCorner(int index) {
    return cornerIterators[index];
  }

  template <class V> hrleSparseOffsetIterator<hrleDomain> &getCorner(V vector) {
    unsigned index = 0;
    for (unsigned i = 0; i < D; ++i) {
      if (vector[i])
        index |= 1 << i;
    }
    return cornerIterators[index];
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  const hrleIndexType &getIndices(unsigned i) { return currentCoords[i]; }

  bool isFinished() const {
    if (compare(currentCoords, maxIndex) > 0) {
      return true;
    } else {
      return false;
    }
  }
};

template <class hrleDomain>
using hrleConstDenseCellIterator = hrleDenseCellIterator<const hrleDomain>;

#endif // HRLE_DENSE_CELL_ITERATOR_HPP
