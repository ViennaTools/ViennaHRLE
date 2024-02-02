#ifndef HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP
#define HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP

#include <cassert>

#include "hrleSparseIterator.hpp"
#include "hrleSparseOffsetIterator.hpp"

/// This class iterates over multiple hrleDomains by using an iterator
/// for each domain and keeping them in sync. A call to next() results
/// in the iterator advancing until the next run in any
/// domain is reached. Therefore, all defined points in multiple
/// domains can be visited sequentially with this iterator. If multiple
/// points are defined at one index (i.e.: the point is defined in more
/// than one domain), it will only stop once.
template <class hrleDomain> class hrleSparseMultiIterator {
public:
  using DomainType = hrleDomain;
  using DomainsType = std::vector<DomainType *>;

private:
  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  static constexpr int D = hrleDomain::dimension;

  DomainsType domains;
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
  hrleSparseMultiIterator<hrleDomain>
  operator++(int) { // use pre increment instead
    return *this;
  }
  hrleSparseMultiIterator<hrleDomain>
  operator--(int) { // use pre decrement instead
    return *this;
  }

public:
  /// A vector with pointers to hrleDomains to iterate over.
  /// The passed hrleVector contains the indices from which to start iterating.
  hrleSparseMultiIterator(DomainsType passedDomains,
                          const hrleVectorType<hrleIndexType, D> v)
      : domains(passedDomains), currentIndices(v) {
    initializeIterators(currentIndices);
  }

  /// A vector with pointers to hrleDomains to iterate over.
  /// The iteration will start from the minimum index of the grid.
  hrleSparseMultiIterator(DomainsType passedDomains)
      : domains(passedDomains),
        currentIndices(passedDomains.back()->getGrid().getMinGridPoint()) {
    initializeIterators(currentIndices);
  }

  /// When this constructor is used, the iteration is started from the minimum
  /// grid index.
  hrleSparseMultiIterator(hrleDomain &passedDomain)
      : currentIndices(passedDomain.getGrid().getMinGridPoint()) {
    domains.push_back(&passedDomain);
    initializeIterators(passedDomain.getGrid().getMinIndex());
  }

  /// Add a domain over which to iterate. This does not affect the current index
  /// in the iteration, so if this should trigger a restart from some start
  /// indices, this has to be done manually.
  void insertNextDomain(hrleDomain &passedDomain) {
    domains.push_back(&passedDomain);
    currentIndices = passedDomain.getGrid().getMinGridPoint();
    iterators.push_back(
        hrleSparseIterator<hrleDomain>(passedDomain, currentIndices));
    goToIndices(currentIndices);
  }

  hrleSparseMultiIterator<hrleDomain> &operator++() {
    next();
    return *this;
  }

  hrleSparseMultiIterator<hrleDomain> &operator--() {
    previous();
    return *this;
  }

  /// go to next defined point in any domain
  void next() {
    if (isFinished()) {
      return;
    }

    // do {
    // find shortest current run to find next run
    hrleVectorType<hrleIndexType, D> endIndices = iterators[0].getEndIndices();
    // std::cout << "End Indices: " << endIndices << std::endl;
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

    // } while (!isDefined() && !isFinished());
  }

  /// got to previous defined point in any domain
  void previous() {
    if (isFinished()) {
      return;
    }

    // do {
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
      it.goToIndices(startIndices);
    }

    currentIndices = startIndices;

    // } while (!isDefined() && !isFinished());
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

  /// Returns a vector of pairs with the domain index and all iterators
  /// currently on defined points.
  std::vector<std::pair<std::size_t, hrleSparseIterator<hrleDomain>>>
  getDefinedIterators() {
    std::vector<std::pair<std::size_t, hrleSparseIterator<hrleDomain>>>
        definedIterators;
    for (std::size_t i = 0; i < iterators.size(); ++i) {
      if (iterators[i].isDefined()) {
        definedIterators.push_back(std::make_pair(i, iterators[i]));
      }
    }
    return definedIterators;
  }

  const DomainsType &getDomains() { return domains; }

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
    currentIndices = v;
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
    currentIndices = v;
  }
};

template <class hrleDomain>
using hrleConstSparseMultiIterator = hrleSparseMultiIterator<const hrleDomain>;

#endif // HRLE_SPARSE_MULTI_DOMAIN_ITERATOR_HPP
