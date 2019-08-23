#ifndef HRLE_CELL_ITERATOR_HPP
#define HRLE_CELL_ITERATOR_HPP

#include <cassert>

#include "hrleOffsetRunsIterator.hpp"

/// This neighbor iterator consists of 2*Dimensions hrleOffsetRunsIterator s
/// for the cartesian neighbors and an hrleRunsIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class hrleCellIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  hrleDomain &domain;
  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleOffsetRunsIterator<hrleDomain>> cornerIterators;

  template <class V> void initialize(const V &v) {
    for (int i = 0; i < (1 << D); ++i) {
      cornerIterators.push_back(hrleOffsetRunsIterator<hrleDomain>(
          domain, BitMaskToVector<D, hrleIndexType>(i), v));
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleCellIterator<hrleDomain> operator++(int); // use pre increment instead
  hrleCellIterator<hrleDomain> operator--(int); // use pre decrement instead

public:
  hrleCellIterator(hrleDomain &passedDomain,
                   const hrleVectorType<hrleIndexType, D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(v);
    if (!isDefined())
      next();
  }

  hrleCellIterator(hrleDomain &passedDomain)
      : domain(passedDomain),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    initialize(passedDomain.getGrid().getMinIndex());
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

  hrleCellIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleCellIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  void next() {
    do {
      const int numCorners = 1 << D;
      std::vector<bool> increment(numCorners, false);
      increment[0] = true;

      hrleVectorType<hrleIndexType, D> end_coords =
          cornerIterators[0].getEndIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (compare(end_coords, cornerIterators[i].getEndIndices())) {
        case 1:
          end_coords = cornerIterators[i].getEndIndices();
          increment = std::vector<bool>(numCorners, false);
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
      std::vector<bool> decrement(numCorners, false);
      decrement[0] = true;

      hrleVectorType<hrleIndexType, D> start_coords =
          cornerIterators[0].getStartIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (compare(start_coords, cornerIterators[i].getStartIndices())) {
        case -1:
          start_coords = cornerIterators[i].getStartIndices();
          decrement = std::vector<bool>(numCorners, false);
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

  hrleOffsetRunsIterator<hrleDomain> &getCorner(unsigned index) {
    return cornerIterators[index];
  }

  hrleOffsetRunsIterator<hrleDomain> &getCorner(int index) {
    return cornerIterators[index];
  }

  template <class V> hrleOffsetRunsIterator<hrleDomain> &getCorner(V vector) {
    unsigned index = 0;
    for (unsigned i = 0; i < D; ++i) {
      if (vector[i])
        index |= 1 << i;
    }
    return cornerIterators[index];
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  const hrleIndexType &getIndices(unsigned i) { return currentCoords[i]; }

  bool isFinished() const { return cornerIterators[0].isFinished(); }
};

template <class hrleDomain>
using hrleConstCellIterator = hrleCellIterator<const hrleDomain>;

#endif // HRLE_CELL_ITERATOR_HPP
