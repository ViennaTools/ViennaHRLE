#ifndef HRLE_DENSE_CELL_ITERATOR_HPP
#define HRLE_DENSE_CELL_ITERATOR_HPP

#include "hrleSparseOffsetIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This neighbor iterator consists of 2*Dimensions SparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class DenseCellIterator {

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  Index<D> currentCoords;
  std::vector<SparseOffsetIterator<hrleDomain>> cornerIterators;
  Index<D> minIndex, maxIndex;

  template <class V> void initialize(const V &v) {
    for (int i = 0; i < (1 << D); ++i) {
      cornerIterators.push_back(SparseOffsetIterator<hrleDomain>(
          domain, BitMaskToVector<D, IndexType>(i), v));
    }
  }

  void incrementIndices(Index<D> &v) {
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

  void decrementIndices(Index<D> &v) {
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
  DenseCellIterator operator++(int) { return *this; }
  // use pre increment instead
  DenseCellIterator operator--(int) { return *this; }
  // use pre decrement instead

public:
  using DomainType = hrleDomain;

  DenseCellIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(currentCoords);
    if (!isDefined())
      next();
  }

  explicit DenseCellIterator(hrleDomain &passedDomain, bool reverse = false)
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

  DenseCellIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  DenseCellIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    const int numCorners = 1 << D;
    // std::vector<bool> increment(numCorners, false);
    // increment[0] = true;

    Index<D> end_coords = currentCoords;
    // cornerIterators[0].getEndIndices();
    for (int i = 0; i < numCorners; i++) {
      switch (Compare(end_coords, cornerIterators[i].getEndIndices())) {
      case 1:
        // TODO
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
    // std::vector<bool> decrement(numCorners, false);
    // decrement[0] = true;

    Index<D> start_coords = currentCoords;
    // cornerIterators[0].getStartIndices();
    for (int i = 0; i < numCorners; i++) {
      switch (Compare(start_coords, cornerIterators[i].getStartIndices())) {
      case -1:
        // TODO
        // start_coords = cornerIterators[i].getStartIndices();
        // decrement = std::vector<bool>(numCorners, false);
      case 0:
        cornerIterators[i].previous();
      }
    }

    decrementIndices(currentCoords);
  }

  SparseOffsetIterator<hrleDomain> &getCorner(unsigned index) {
    return cornerIterators[index];
  }

  SparseOffsetIterator<hrleDomain> &getCorner(int index) {
    return cornerIterators[index];
  }

  template <class V> SparseOffsetIterator<hrleDomain> &getCorner(V vector) {
    unsigned index = 0;
    for (unsigned i = 0; i < D; ++i) {
      if (vector[i])
        index |= 1 << i;
    }
    return cornerIterators[index];
  }

  const Index<D> &getIndices() { return currentCoords; }

  const IndexType &getIndices(unsigned i) { return currentCoords[i]; }

  const DomainType &getDomain() { return domain; }

  bool isFinished() const {
    if (compare(currentCoords, maxIndex) > 0) {
      return true;
    } else {
      return false;
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
        ++(*this);
      }
    } else {
      while (v < currentCoords) {
        --(*this);
      }
    }
  }
};

template <class hrleDomain>
using ConstDenseCellIterator = DenseCellIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_DENSE_CELL_ITERATOR_HPP
