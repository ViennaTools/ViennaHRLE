#ifndef HRLE_CELL_ITERATOR_HPP
#define HRLE_CELL_ITERATOR_HPP

#include <array>
#include <utility>

#include "hrleSparseOffsetIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This neighbor iterator consists of 2*Dimensions SparseOffsetIterator s
/// for the cartesian neighbors and an hrleSparseIterator
/// for the center.
/// Whenever one of these (2*Dimensions+1) iterators reach a defined grid point,
/// the iterator stops.
template <class hrleDomain> class SparseCellIterator {
public:
  using DomainType = hrleDomain;
  using OffsetIterator = SparseOffsetIterator<hrleDomain>;

private:
  using ValueType = std::conditional_t<std::is_const_v<hrleDomain>,
                                       const typename hrleDomain::ValueType,
                                       typename hrleDomain::ValueType>;
  static constexpr int D = hrleDomain::dimension;
  static constexpr int numCorners = 1 << D;

  hrleDomain &domain;
  Index<D> currentCoords;
  std::array<OffsetIterator, numCorners> cornerIterators;

  template <std::size_t... Is>
  static std::array<OffsetIterator, numCorners>
  makeCornerIteratorsImpl(hrleDomain &passedDomain, const Index<D> &v,
                          std::index_sequence<Is...>) {
    return {OffsetIterator(passedDomain,
                           BitMaskToIndex<D>(static_cast<unsigned>(Is)), v)...};
  }

  static std::array<OffsetIterator, numCorners>
  makeCornerIterators(hrleDomain &passedDomain, const Index<D> &v) {
    return makeCornerIteratorsImpl(
        passedDomain, v,
        std::make_index_sequence<static_cast<std::size_t>(numCorners)>{});
  }

public:
  explicit SparseCellIterator(hrleDomain &passedDomain)
      : domain(passedDomain), currentCoords(domain.getGrid().getMinGridPoint()),
        cornerIterators(makeCornerIterators(
            passedDomain, passedDomain.getGrid().getMinIndex())) {
    if (!isDefined())
      next();
  }

  SparseCellIterator(hrleDomain &passedDomain, const Index<D> &v)
      : domain(passedDomain), currentCoords(v),
        cornerIterators(makeCornerIterators(passedDomain, v)) {
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
    for (int i = 0; i < numCorners; i++) {
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
    static_assert(numCorners <= 8);

    do {
      std::uint_fast8_t incrementMask = 1u; // corner 0

      auto end_coords = cornerIterators[0].getEndIndices();

      for (int i = 1; i < numCorners; ++i) {
        const auto &cornerEnd = cornerIterators[i].getEndIndices();
        const int cmp = Compare(end_coords, cornerEnd);

        if (cmp > 0) {
          end_coords = cornerEnd;
          incrementMask = static_cast<std::uint_fast8_t>(1u << i);
        } else if (cmp == 0) {
          incrementMask |= static_cast<std::uint_fast8_t>(1u << i);
        }
      }

      for (int i = 0; i < numCorners; ++i) {
        if (incrementMask & (1u << i))
          cornerIterators[i].next();
      }

      currentCoords = domain.getGrid().incrementIndices(end_coords);
    } while (!isDefined() && !isFinished());
  }

  void previous() {
    static_assert(numCorners <= 8);

    do {
      std::uint_fast8_t decrementMask = 1u; // corner 0

      auto start_coords = cornerIterators[0].getStartIndices();

      for (int i = 1; i < numCorners; ++i) {
        const auto &cornerStart = cornerIterators[i].getStartIndices();
        const int cmp = Compare(start_coords, cornerStart);

        if (cmp < 0) {
          start_coords = cornerStart;
          decrementMask = static_cast<std::uint_fast8_t>(1u << i);
        } else if (cmp == 0) {
          decrementMask |= static_cast<std::uint_fast8_t>(1u << i);
        }
      }

      for (int i = 0; i < numCorners; ++i) {
        if (decrementMask & (1u << i))
          cornerIterators[i].previous();
      }

      currentCoords = domain.getGrid().decrementIndices(start_coords);
    } while (!isDefined() && !isFinished());
  }

  OffsetIterator &getCorner(unsigned index) { return cornerIterators[index]; }

  OffsetIterator const &getCorner(unsigned index) const {
    return cornerIterators[index];
  }

  OffsetIterator &getCorner(int index) { return cornerIterators[index]; }

  OffsetIterator const &getCorner(int index) const {
    return cornerIterators[index];
  }

  OffsetIterator &getCorner(const Index<D> &vector) {
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

  bool isFinished() const { return cornerIterators[0].isFinished(); }
};

template <class hrleDomain>
using ConstSparseCellIterator = SparseCellIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_CELL_ITERATOR_HPP
