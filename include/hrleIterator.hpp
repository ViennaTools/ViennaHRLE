#ifndef HRLE_ITERATOR_HPP
#define HRLE_ITERATOR_HPP

#include "hrleRunsIterator.hpp"

/// This iterator iterates over the entire structure and stops at every single
/// grid point. It also stops on undefined grid points. Therefore, it can be
/// used to easily generate a full grid from the sparse data set.
template <class hrleDomain> class hrleIterator {

  static constexpr int D = hrleDomain::dimension;
  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  hrleDomain &domain;
  hrleRunsIterator<hrleDomain> runsIterator;
  hrleVectorType<hrleIndexType, D> currentIndices;
  hrleVectorType<hrleIndexType, D> minIndex;
  hrleVectorType<hrleIndexType, D> maxIndex;

  void incrementIndices(hrleVectorType<hrleIndexType, D> &v) {
    int dim = 0;
    for (; dim < D - 1; ++dim) {
      bool posInfinite = domain.getGrid().isPosBoundaryInfinite(dim);
      bool negInfinite = domain.getGrid().isNegBoundaryInfinite(dim);
      if (v[dim] < (posInfinite ? domain.getMaxRunBreak(dim)
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
      v[dim] = (posInfinite ? domain.getMaxRunBreak(dim)
                            : domain.getGrid().getMaxGridPoint(dim));
    }
    --v[dim];
  }

public:
  hrleIterator(hrleDomain &passedDomain, bool reverse = false)
      : domain(passedDomain), runsIterator(passedDomain, reverse) {
    auto &grid = domain.getGrid();
    for(unsigned i=0; i<D; ++i){
      minIndex[i] = (grid.isNegBoundaryInfinite(i))?domain.getMinRunBreak(i):grid.getMinBounds(i);
      maxIndex[i] = (grid.isPosBoundaryInfinite(i))?domain.getMaxRunBreak(i):grid.getMaxBounds(i);
    }
    if (reverse)
      currentIndices = maxIndex;
    else
      currentIndices = minIndex;
  }

  template <class V>
  hrleIterator(hrleDomain &passedDomain, V &v)
      : domain(passedDomain), runsIterator(passedDomain, v) {}

  hrleIterator<hrleDomain> &operator++() {
    // move iterator with currentIndices if they are the same
    switch (compare(runsIterator.getEndIndices(), currentIndices)) {
    case -1:
    case 0:
      if (!runsIterator.isFinished())
        runsIterator.next();
    default:
      incrementIndices(currentIndices);
    }
    return *this;
  }

  hrleIterator<hrleDomain> operator++(int) {
    hrleIterator<hrleDomain> tmp(*this);
    ++(*this);
    return tmp;
  }

  hrleIterator<hrleDomain> &operator--() {
    // move iterator with currentIndices if they are the same
    switch (compare(runsIterator.getStartIndices(), currentIndices)) {
    case 1:
    case 0:
      if (!runsIterator.isFinished())
        runsIterator.previous();
    default:
      decrementIndices(currentIndices);
    }
    return *this;
  }

  hrleIterator<hrleDomain> operator--(int) {
    hrleIterator<hrleDomain> tmp(*this);
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

  void print() {
    std::cout << currentIndices << std::endl;
    runsIterator.print();
  }
};

template <class hrleDomain>
using hrleConstIterator = hrleIterator<const hrleDomain>;

#endif // HRLE_ITERATOR_HPP
