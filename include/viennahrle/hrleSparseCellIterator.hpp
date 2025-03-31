#ifndef HRLE_CELL_ITERATOR_HPP
#define HRLE_CELL_ITERATOR_HPP

#include "hrleSparseOffsetIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This neighbor iterator consists of 2*Dimensions SparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class SparseCellIterator {

  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  static constexpr int D = hrleDomain::dimension;
  static constexpr int numCorners = 1 << D;

  hrleDomain &domain;
  Index<D> currentCoords;
  std::vector<SparseOffsetIterator<hrleDomain>> cornerIterators;

  template <class V> void initialize(const V &v) {
    for (unsigned i = 0; i < 1 << D; ++i) {
      cornerIterators.push_back(
          SparseOffsetIterator<hrleDomain>(domain, BitMaskToIndex<D>(i), v));
    }
  }

public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

  SparseCellIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v) {

    initialize(v);
    if (!isDefined())
      next();
  }

  explicit SparseCellIterator(hrleDomain &passedDomain)
      : domain(passedDomain),
        currentCoords(domain.getGrid().getMinGridPoint()) {

    initialize(passedDomain.getGrid().getMinIndex());
    if (!isDefined())
      next();
  }

  // delete post in/decrement, since they should not be used, due to the
  // size of the structure
  SparseCellIterator operator++(int) = delete; // use pre increment instead
  SparseCellIterator operator--(int) = delete; // use pre decrement instead

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

  SparseCellIterator &operator++() {
    next();
    return *this;
  }

  SparseCellIterator &operator--() {
    previous();
    return *this;
  }

  void next() {
    do {
      std::array<bool, numCorners> increment;
      increment.fill(false);
      increment[0] = true;

      auto end_coords = cornerIterators[0].getEndIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (Compare(end_coords, cornerIterators[i].getEndIndices())) {
        case 1:
          end_coords = cornerIterators[i].getEndIndices();
          increment.fill(false);
        case 0:
          increment[i] = true;
        default:
          break;
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
      std::array<bool, numCorners> decrement;
      decrement.fill(false);
      decrement[0] = true;

      auto start_coords = cornerIterators[0].getStartIndices();
      for (int i = 1; i < numCorners; i++) {
        switch (Compare(start_coords, cornerIterators[i].getStartIndices())) {
        case -1:
          start_coords = cornerIterators[i].getStartIndices();
          decrement.fill(false);
        case 0:
          decrement[i] = true;
        default:
          break;
        }
      }

      for (int i = 0; i < numCorners; ++i) {
        if (decrement[i])
          cornerIterators[i].previous();
      }
      currentCoords = domain.getGrid().decrementIndices(start_coords);
    } while (!isDefined() && !isFinished());
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
    assert(index < numCorners);
    return cornerIterators[index];
  }

  const Index<D> &getIndices() { return currentCoords; }

  const IndexType &getIndices(unsigned i) { return currentCoords[i]; }

  const DomainType &getDomain() { return domain; }

  bool isFinished() const { return cornerIterators[0].isFinished(); }
};

template <class hrleDomain>
using ConstSparseCellIterator = SparseCellIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_CELL_ITERATOR_HPP
