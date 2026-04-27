#ifndef HRLE_DENSE_CELL_ITERATOR_HPP
#define HRLE_DENSE_CELL_ITERATOR_HPP

#include "hrleSparseOffsetIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This iterator consists of 2^Dimensions SparseOffsetIterators
/// for the cartesian neighbors and an SparseIterator for the center.
template <class hrleDomain> class DenseCellIterator {
public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

private:
  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  static constexpr int D = hrleDomain::dimension;
  static constexpr int numCorners = 1 << D;

  hrleDomain &domain;
  Index<D> currentCoords;
  std::vector<OffsetIterator> cornerIterators;
  Index<D> minIndex, maxIndex;

  template <class V> void initialize(const V &v) {
    cornerIterators.reserve(numCorners);
    for (unsigned i = 0; i < numCorners; ++i) {
      cornerIterators.emplace_back(domain, BitMaskToIndex<D>(i), v);
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

public:
  DenseCellIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(currentCoords);
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
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  DenseCellIterator operator++(int) = delete; // use pre increment instead
  DenseCellIterator operator--(int) = delete; // use pre decrement instead

  bool isDefined() const {
    for (unsigned i = 0; i < D; ++i) {
      if (!domain.getGrid().isBoundaryPeriodic(i) &&
          currentCoords[i] == domain.getGrid().getMaxGridPoint(i)) {
        return false;
      }
    }
    for (int i = 0; i < numCorners; i++) {
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
      default:
        break;
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
      default:
        break;
      }
    }

    decrementIndices(currentCoords);
  }

  OffsetIterator &getCorner(unsigned index) { return cornerIterators[index]; }

  OffsetIterator const &getCorner(unsigned index) const {
    return cornerIterators[index];
  }

  OffsetIterator &getCorner(int index) { return cornerIterators[index]; }

  OffsetIterator const &getCorner(int index) const {
    return cornerIterators[index];
  }

  template <class V> OffsetIterator &getCorner(V vector) {
    unsigned index = 0;
    for (unsigned i = 0; i < D; ++i) {
      if (vector[i])
        index |= 1 << i;
    }
    assert(index < numCorners);
    return cornerIterators[index];
  }

  const Index<D> &getIndices() const { return currentCoords; }

  const IndexType &getIndices(unsigned i) const { return currentCoords[i]; }

  const DomainType &getDomain() const { return domain; }

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
