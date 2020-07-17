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

  std::vector<hrleDomain *> domains;
  hrleVectorType<hrleIndexType, D> currentIndices;
  std::vector<hrleSparseIterator<hrleDomain>> iterators;

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
  /// A vector with pointers to hrleDomains to iterate over.
  /// The passed hrleVector contains the indices from which to start iterating.
  hrleSparseMultiDomainIterator(std::vector<hrleDomain *> passedDomains,
                                const hrleVectorType<hrleIndexType, D> v)
      : domains(passedDomains), currentIndices(v) {
    initializeIterators(currentIndices);
  }

  /// A vector with pointers to hrleDomains to iterate over.
  /// The iteration will start from the minimum index of the grid.
  hrleSparseMultiDomainIterator(std::vector<hrleDomain *> passedDomains)
      : domains(passedDomains),
        currentIndices(*(passedDomains[0]).getGrid().getGridMin()) {
    initializeIterators(currentIndices);
  }

  /// When this constructor is used, the iteration is started from the minimum
  /// grid index.
  hrleSparseMultiDomainIterator(hrleDomain &passedDomain)
      : currentIndices(passedDomain.getGrid().getMinGridPoint()) {
    domains.push_back(&passedDomain);
    initializeIterators(passedDomain.getGrid().getMinIndex());
  }

  /// Add a domain over which to iterate. This does not affect the current index
  /// in the iteration, so if this should trigger a restart from some start
  /// indices, this has to be done manually.
  void insertNextDomain(hrleDomain &passedDomain) {
    domains.push_back(&passedDomain);
    iterators.push_back(
        hrleSparseIterator<hrleDomain>(passedDomain, currentIndices));
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
    if (isFinished()) {
      return;
    }

    do {
      // find shortest current run to find next run
      hrleVectorType<hrleIndexType, D> endIndices =
          iterators[0].getEndIndices();
      for (unsigned i = 1; i < iterators.size(); ++i) {
        if (!iterators[i].isFinished() &&
            iterators[i].getEndIndices() < endIndices) {
          endIndices = iterators[i].getEndIndices();
        }
      }

      endIndices = domains[0]->getGrid().incrementIndices(endIndices);

      // now advance all iterators to reach next defined run
      for (auto &it : iterators) {
        it.goToIndicesSequential(endIndices);
      }

      currentIndices = endIndices;

    } while (!isDefined() && !isFinished());
  }

  /// got to previous defined point in any domain
  void previous() {
    if (isFinished()) {
      return;
    }

    do {
      // find shortest current run to find next run
      hrleVectorType<hrleIndexType, D> startIndices =
          iterators[0].getStartIndices();
      for (unsigned i = 1; i < iterators.size(); ++i) {
        if (!iterators[i].isFinished() &&
            iterators[i].getStartIndices() < startIndices) {
          startIndices = iterators[i].getStartIndices();
        }
      }

      startIndices = domains[0]->getGrid().decrementIndices(startIndices);

      // now advance all iterators to reach next defined run
      for (auto &it : iterators) {
        it.goToIndicesSequential(startIndices);
      }

      currentIndices = startIndices;

    } while (!isDefined() && !isFinished());
  }

  /// Returns the iterator used for the domain at index. This is the index
  /// the domain had in the std::vector upon initialisation, or the index
  /// it obtained through calls to insertNextDomain.
  hrleSparseIterator<hrleDomain> &getIterator(int index) {
    return iterators[index];
  }

  hrleIndexType getIndex(int i) { return currentIndices[i]; }

  const hrleVectorType<hrleIndexType, D> &getIndices() {
    return currentIndices;
  }

  std::size_t getNumberOfDomains() { return domains.size(); }

  /// Returns whether there is any defined point at the current index.
  bool isDefined() const {
    for (unsigned i = 0; i < iterators.size(); i++) {
      if (iterators[i].isDefined())
        return true;
    }
    return false;
  }

  /// Returns a vector with all iterators currently on defined points.
  std::vector<hrleSparseIterator<hrleDomain>> getDefinedIterators() {
    std::vector<hrleSparseIterator<hrleDomain>> definedIterators;
    for (auto &it : iterators) {
      if (it.isDefined()) {
        definedIterators.push_back(it);
      }
    }
    return definedIterators;
  }

  /// If all iterators in all domains are finished, this will return true.
  bool isFinished() const {
    for (auto &it : iterators) {
      if (!it.isFinished()) {
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
    if (v >= currentIndices) {
      while (v > currentIndices) {
        next();
      }
    } else {
      while (v < currentIndices) {
        previous();
      }
    }
  }
};

template <class hrleDomain>
using hrleConstSparseMultiDomainIterator =
    hrleSparseMultiDomainIterator<const hrleDomain>;

#endif // HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP
