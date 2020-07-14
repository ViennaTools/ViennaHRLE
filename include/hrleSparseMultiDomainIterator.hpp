#ifndef HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP
#define HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP

#include <cassert>

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"

/// This class iterates over multiple hrleDomains by using an iterator
/// for each domain and keeping them in sync. A call to next() results
/// in the iterator advancing until the next defined point in any
/// domain is reached. Therefore, all defined points in multiple
/// domains can be visited sequentially with this iterator. If multiple
/// points are defined at one index (i.e.: the point is defined in more
/// than one domain), it will only stop once.
template <class hrleDomain> class hrleSparseMultiDomainIterator {

  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  std::vector<hrleDomain*> domains;
  hrleVectorType<hrleIndexType, D> currentCoords;
  std::vector<hrleSparseIterator<hrleDomain>> iterators;

  bool isDefined() const {
    for (unsigned i = 0; i < iterators.size(); i++) {
      if (iterators[i].isDefined())
        return true;
    }
    return false;
  }

  template <class V> void initializeIterators(const V &v) {
    iterators.clear();
    for (unsigned i = 0; i < domains.size(); ++i) {
      iterators.push_back(hrleSparseIterator<hrleDomain>(*(domains[i]), v));
    }
  }

  // make post in/decrement private, since they should not be used, due to the
  // size of the structure
  hrleSparseMultiDomainIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleSparseMultiDomainIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  hrleSparseMultiDomainIterator(std::vector<hrleDomain*> passedDomains,
                         const hrleVectorType<hrleIndexType, D> &v)
      : domains(passedDomains), currentCoords(v) {
    initializeIterators(v);
  }

  hrleSparseMultiDomainIterator(hrleDomain &passedDomain)
      : currentCoords(passedDomain.getGrid().getMinGridPoint()) {
    domains.push_back(&passedDomain);
    initializeIterators(passedDomain.getGrid().getMinIndex());
  }

  hrleSparseMultiDomainIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleSparseMultiDomainIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  /// go to next defined point in any domain
  void next() {
    if(isFinished) {
      return;
    }

    do {
      // find shortest current run to find next run
      hrleVectorType<hrleIndexType, D> end_coords = iterators[0].getEndIndices();
      for(unsigned i = 1; i < iterators.size(); ++i) {
        if(iterators[i].getEndIndices() < end_coords) {
          end_coords = iterators[i].getEndIndices();
        }
      }

      // now advance all iterators to reach next defined run


    } while(!isDefined() && !isFinished());
    


    // const int numNeighbors = 2 * order * D;
    // std::vector<bool> increment(numNeighbors + 1, false);
    // increment[numNeighbors] = true;

    // hrleVectorType<hrleIndexType, D> end_coords =
    //     centerIterator.getEndIndices();
    // for (int i = 0; i < numNeighbors; i++) {
    //   switch (compare(end_coords, neighborIterators[i].getEndIndices())) {
    //   case 1:
    //     end_coords = neighborIterators[i].getEndIndices();
    //     increment = std::vector<bool>(numNeighbors + 1, false);
    //   case 0:
    //     increment[i] = true;
    //   }
    // }

    // if (increment[numNeighbors])
    //   centerIterator.next();
    // for (int i = 0; i < numNeighbors; ++i)
    //   if (increment[i])
    //     neighborIterators[i].next();

    // currentCoords = domain.getGrid().incrementIndices(end_coords);
  }

  void previous() {
    // const int numNeighbors = 2 * order * D;
    // std::vector<bool> decrement(numNeighbors + 1, false);
    // decrement[numNeighbors] = true;

    // hrleVectorType<hrleIndexType, D> start_coords =
    //     centerIterator.getStartIndices();
    // for (int i = 0; i < numNeighbors; i++) {
    //   switch (compare(start_coords, neighborIterators[i].getStartIndices())) {
    //   case -1:
    //     start_coords = neighborIterators[i].getStartIndices();
    //     decrement = std::vector<bool>(numNeighbors + 1, false);
    //   case 0:
    //     decrement[i] = true;
    //   }
    // }

    // if (decrement[numNeighbors])
    //   centerIterator.previous();
    // for (int i = 0; i < numNeighbors; ++i) {
    //   if (decrement[i])
    //     neighborIterators[i].previous();
    // }
    // currentCoords = domain.getGrid().decrementIndices(start_coords);
  }

  /// Returns the iterator used for the domain at index. This is the index
  /// the domain had in the std::vector upon initialisation, or the index
  /// it obtained through calls to insertNextDomain.
  hrleSparseOffsetIterator<hrleDomain> &getDomainIterator(int index) {
    return iterators[index];
  }

  const hrleVectorType<hrleIndexType, D> &getIndices() { return currentCoords; }

  bool isFinished() const {
    for(auto &it : iterators) {
      if(!it.isFinished()) {
        return false;
      }
    }
    return true;
  }

  /// Sets the iterator to position v.
  /// Uses random access to move, so it is be slower
  /// than goToIndicesSequential for repeated serial calls.
  template <class V> void goToIndices(V &v) {
    for (auto &it : iterators) {
      it.goToIndices(v);
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
        next();
      }
    } else {
      while (v < currentCoords) {
        previous();
      }
    }
  }
};

template <class hrleDomain>
using hrleConstSparseStarIterator = hrleSparseMultiDomainIterator<const hrleDomain>;

#endif // HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP
