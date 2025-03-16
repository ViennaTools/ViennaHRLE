#ifndef HRLE_ITERATOR_HPP
#define HRLE_ITERATOR_HPP

#include "hrleSparseIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This iterator iterates over the entire structure and stops at every single
/// grid point. It also stops on undefined grid points. Therefore, it can be
/// used to easily generate a full grid from the sparse data set.
template <class hrleDomain> class DenseIterator {

  static constexpr int D = hrleDomain::dimension;
  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  hrleDomain &domain;
  SparseIterator<hrleDomain> runsIterator;
  Index<D> currentIndices;
  Index<D> minIndex;
  Index<D> maxIndex;

  void incrementIndices(Index<D> &v) {
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

  void decrementIndices(Index<D> &v) {
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

  explicit DenseIterator(hrleDomain &passedDomain, bool reverse = false)
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
  DenseIterator(hrleDomain &passedDomain, V &v)
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

  DenseIterator<hrleDomain> &operator++() {
    // move iterator with currentIndices if they are the same
    incrementIndices(currentIndices);
    while (runsIterator.getEndIndices() < currentIndices &&
           !runsIterator.isFinished()) {
      runsIterator.next();
    }

    return *this;
  }

  DenseIterator<hrleDomain> operator++(int) {
    DenseIterator<hrleDomain> tmp(*this);
    ++(*this);
    return tmp;
  }

  DenseIterator<hrleDomain> &operator--() {
    // move iterator with currentIndices if they are the same
    decrementIndices(currentIndices);
    while (runsIterator.getEndIndices() > currentIndices &&
           !runsIterator.isFinished()) {
      runsIterator.previous();
    }

    return *this;
  }

  DenseIterator<hrleDomain> operator--(int) {
    DenseIterator<hrleDomain> tmp(*this);
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
    if (Compare(currentIndices, minIndex) < 0) {
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
    if (Compare(currentIndices, maxIndex) > 0)
      return true;

    return false;
  }

  bool isDefined() const { return runsIterator.isDefined(); }

  SizeType getPointId() const { return runsIterator.getPointId(); }

  ValueType &getValue() { return runsIterator.getValue(); }

  IndexType getIndex(int dimension) { return currentIndices[dimension]; }

  Index<D> getIndices() { return currentIndices; }

  Index<D> getIteratorIndices() { return runsIterator.getStartIndices(); }

  const DomainType &getDomain() { return domain; }

  void print() {
    std::cout << currentIndices << std::endl;
    runsIterator.print();
  }
};

template <class hrleDomain>
using ConstDenseIterator = DenseIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_ITERATOR_HPP
