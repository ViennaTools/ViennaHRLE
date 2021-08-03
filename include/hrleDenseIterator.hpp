#ifndef HRLE_ITERATOR_HPP
#define HRLE_ITERATOR_HPP

#include "hrleSparseIterator.hpp"

/// This iterator iterates over the entire structure and stops at every single
/// grid point. It also stops on undefined grid points. Therefore, it can be
/// used to easily generate a full grid from the sparse data set.
template <class hrleDomain> class hrleDenseIterator {

  static constexpr int D = hrleDomain::dimension;
  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  hrleDomain &domain;
  hrleSparseIterator<hrleDomain> runsIterator;
  hrleVectorType<hrleIndexType, D> currentIndices;
  hrleVectorType<hrleIndexType, D> minIndex;
  hrleVectorType<hrleIndexType, D> maxIndex;

  void incrementIndices(hrleVectorType<hrleIndexType, D> &v) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      bool posInfinite = domain.getGrid().isPosBoundaryInfinite(dim);
      bool negInfinite = domain.getGrid().isNegBoundaryInfinite(dim);
      if (v[dim] < (posInfinite ? domain.getMaxRunBreak(dim) - 1
                                : domain.getGrid().getMaxGridPoint(dim)))
        break;
      v[dim] = (negInfinite ? domain.getMinRunBreak(dim)
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
      v[dim] = (posInfinite ? domain.getMaxRunBreak(dim) - 1
                            : domain.getGrid().getMaxGridPoint(dim));
    }
    --v[dim];
  }

public:
  using DomainType = hrleDomain;

  hrleDenseIterator(hrleDomain &passedDomain, bool reverse = false)
      : domain(passedDomain), runsIterator(passedDomain, reverse) {
    auto &grid = domain.getGrid();
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = (grid.isNegBoundaryInfinite(i)) ? domain.getMinRunBreak(i)
                                                    : grid.getMinBounds(i);
      maxIndex[i] = (grid.isPosBoundaryInfinite(i))
                        ? domain.getMaxRunBreak(i) - 1
                        : grid.getMaxBounds(i);
    }
    if (reverse)
      currentIndices = maxIndex;
    else
      currentIndices = minIndex;
    runsIterator.goToIndicesSequential(currentIndices);
  }

  template <class V>
  hrleDenseIterator(hrleDomain &passedDomain, V &v)
      : domain(passedDomain), runsIterator(passedDomain, v) {
    auto &grid = domain.getGrid();
    for (unsigned i = 0; i < D; ++i) {
      minIndex[i] = (grid.isNegBoundaryInfinite(i)) ? domain.getMinRunBreak(i)
                                                    : grid.getMinBounds(i);
      maxIndex[i] = (grid.isPosBoundaryInfinite(i))
                        ? domain.getMaxRunBreak(i) - 1
                        : grid.getMaxBounds(i);
    }
    currentIndices = minIndex;
  }

  hrleDenseIterator<hrleDomain> &operator++() {
    // move iterator with currentIndices if they are the same
    incrementIndices(currentIndices);
    while (runsIterator.getEndIndices() < currentIndices &&
           !runsIterator.isFinished()) {
      runsIterator.next();
    }

    return *this;
  }

  hrleDenseIterator<hrleDomain> operator++(int) {
    hrleDenseIterator<hrleDomain> tmp(*this);
    ++(*this);
    return tmp;
  }

  hrleDenseIterator<hrleDomain> &operator--() {
    // move iterator with currentIndices if they are the same
    decrementIndices(currentIndices);
    while (runsIterator.getEndIndices() > currentIndices &&
           !runsIterator.isFinished()) {
      runsIterator.previous();
    }

    return *this;
  }

  hrleDenseIterator<hrleDomain> operator--(int) {
    hrleDenseIterator<hrleDomain> tmp(*this);
    --(*this);
    return tmp;
  }

  // safe version of operator++, returns false if iterator is already done
  bool next() {
    // if max index is reached, iterator is done
    if (isFinished())
      return false;
    ++(*this);
    return true;
  }

  // safe version of operator--, returns false if iterator is already at the
  // start
  bool previous() {
    // if min index is reached, iterator is done
    if (compare(minIndex) < 0) {
      return false;
    }
    ++(*this);
    return true;
  }

  template <class V> void goToIndices(V &v) {
    currentIndices = v;
    runsIterator.goToIndices(v);
  }

  bool isFinished() {
    if (compare(currentIndices, maxIndex) > 0) {
      return true;
    } else {
      return false;
    }
  }

  hrleValueType &getValue() { return runsIterator.getValue(); }

  hrleIndexType getIndex(int dimension) { return currentIndices[dimension]; }

  hrleVectorType<hrleIndexType, D> getIndices() { return currentIndices; }

  hrleVectorType<hrleIndexType, D> getIteratorIndices() {
    return runsIterator.getStartIndices();
  }

  const DomainType &getDomain() { return domain; }

  void print() {
    std::cout << currentIndices << std::endl;
    runsIterator.print();
  }
};

template <class hrleDomain>
using hrleConstDenseIterator = hrleDenseIterator<const hrleDomain>;

#endif // HRLE_ITERATOR_HPP
