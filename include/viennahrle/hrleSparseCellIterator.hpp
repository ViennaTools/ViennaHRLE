#ifndef HRLE_CELL_ITERATOR_HPP
#define HRLE_CELL_ITERATOR_HPP

#include "hrleSparseOffsetIterator.hpp"

/// This neighbor iterator consists of 2*Dimensions hrleSparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class hrleSparseCellIterator {

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::hrleValueType,
                             typename hrleDomain::hrleValueType>
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleSparseOffsetIterator<hrleDomain>> cornerIterators;

  template <class V> void initialize(const V &v) {
    for (int i = 0; i < (1 << D); ++i) {
      cornerIterators.push_back(hrleSparseOffsetIterator<hrleDomain>(
          domain, hrleUtil::BitMaskToVector<D, hrleIndexType>(i), v));
    }
  }

public:
  using DomainType = hrleDomain;

  hrleSparseCellIterator(hrleDomain &passedDomain,
                         const hrleVectorType<hrleIndexType, D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(v);
    if (!isDefined())
      next();
  }

  explicit hrleSparseCellIterator(hrleDomain &passedDomain)
      : domain(passedDomain),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    initialize(passedDomain.getGrid().getMinIndex());
    if (!isDefined())
      next();
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  hrleSparseCellIterator operator++(int) = delete; // use pre increment instead
  hrleSparseCellIterator operator--(int) = delete; // use pre decrement instead

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

  hrleSparseCellIterator &operator++() {
    next();
    return *this;
  }

  hrleSparseCellIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    do {
      const int numCorners = 1 << D;
      std::array<bool, 1 << D> increment;
      increment.fill(false);
      increment[0] = true;

      hrleVectorType<hrleIndexType, D> end_coords =
          cornerIterators[0].getEndIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (
            hrleUtil::Compare(end_coords, cornerIterators[i].getEndIndices())) {
        case 1:
          end_coords = cornerIterators[i].getEndIndices();
          increment.fill(false);
        case 0:
          increment[i] = true;
        }
      }

      for (int i = 0; i < numCorners; ++i)
        if (increment[i])
          cornerIterators[i].next();

      currentCoords = domain.getGrid().incrementIndices(end_coords);
    } while (!isDefined() && !isFinished());
  }

  void previous() {
    do {
      const int numCorners = 1 << D;
      std::array<bool, 1 << D> decrement;
      decrement.fill(false);
      decrement[0] = true;

      hrleVectorType<hrleIndexType, D> start_coords =
          cornerIterators[0].getStartIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (hrleUtil::Compare(start_coords,
                                  cornerIterators[i].getStartIndices())) {
        case -1:
          start_coords = cornerIterators[i].getStartIndices();
          decrement.fill(false);
        case 0:
          decrement[i] = true;
        }
      }

      for (int i = 0; i < numCorners; ++i) {
        if (decrement[i])
          cornerIterators[i].previous();
      }
      currentCoords = domain.getGrid().decrementIndices(start_coords);
    } while (!isDefined() && !isFinished());
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

  const DomainType &getDomain() { return domain; }

  bool isFinished() const { return cornerIterators[0].isFinished(); }
};

template <class hrleDomain>
using hrleConstSparseCellIterator = hrleSparseCellIterator<const hrleDomain>;

#endif // HRLE_CELL_ITERATOR_HPP
