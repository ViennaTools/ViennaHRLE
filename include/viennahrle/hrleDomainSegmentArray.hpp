#ifndef HRLE_DOMAIN_SEGMENT_ARRAY_HPP
#define HRLE_DOMAIN_SEGMENT_ARRAY_HPP

#include "hrleDomainSegment.hpp"
#include "hrleGrid.hpp"

namespace viennahrle {
template <class T = double, int D = 3> class DomainSegmentArray {
  // contains the vector of hrleDomainSegment items that are used for
  // parallelization

  typedef DomainSegment<T, D> *segmentPointer; // pointer to memory containing
                                               // hrleDomainSegment object
  std::vector<segmentPointer> domainSegmentPointers;

public:
  // assignment operator and copy-constructor are deleted since vector only
  // holds pointers, use deepCopy or shallowCopy instead
  DomainSegmentArray &operator=(const DomainSegmentArray &s) = delete;
  DomainSegmentArray(const DomainSegmentArray &s) = delete;

  using size_type = typename std::vector<segmentPointer>::size_type;

  const DomainSegment<T, D> &operator[](size_type i) const {
    return *(domainSegmentPointers[i]);
  }

  DomainSegment<T, D> &operator[](size_type i) {
    return *(domainSegmentPointers[i]);
  }

  const DomainSegment<T, D> &back() const {
    return *(domainSegmentPointers.back());
  }

  DomainSegment<T, D> &back() { return *(domainSegmentPointers.back()); }

  void swap(DomainSegmentArray &s) noexcept {
    domainSegmentPointers.swap(s.domainSegmentPointers);
  }

  void clear() {
    for (auto &segP : domainSegmentPointers)
      if (segP)
        delete segP;
    domainSegmentPointers.clear();
  }

  /// initialize the parallelization by implementing i number of
  /// hrleDomainSegment objects through a segmentPointer object, that
  /// each correspond to a parallel core
  void
  initialize(size_type i, const Grid<D> &g,
             AllocationType<SizeType, D> a = AllocationType<SizeType, D>()) {
    clear(); // delete all data in the array

    domainSegmentPointers.resize(i);

    a *= a.allocationFactor; // increase allocation by extra margin defined in
                             // AllocationType
    a /= i;

#pragma omp parallel for schedule(static, 1)
    // parallelization - Iterations divided into chunks of size 1.
    // Each chunk is assigned to a thread
    for (int k = 0; k < static_cast<int>(domainSegmentPointers.size()); ++k) {
      domainSegmentPointers[k] = segmentPointer(new DomainSegment<T, D>(g, a));
    }
  }

  size_type getNumberOfSegments() const {
    // return the size of the vector for parallelization
    return domainSegmentPointers.size();
  }

  void deepCopy(const Grid<D> *g, const DomainSegmentArray &s) {
    clear(); // deallocate all domain segments
    domainSegmentPointers.resize(s.getNumberOfSegments());

#pragma omp parallel for schedule(static, 1)
    // parallelization - Iterations divided into chunks of size 1.
    // Each chunk is assigned to a thread
    for (int k = 0; k < static_cast<int>(s.getNumberOfSegments()); ++k) {
      domainSegmentPointers[k] =
          segmentPointer(new DomainSegment<T, D>(*g, s[k]));
    }
  }

  void shallowCopy(const DomainSegmentArray &s) {
    for (int k = 0; k < static_cast<int>(s.getNumberOfSegments()); ++k) {
      domainSegmentPointers[k] = segmentPointer(&s[k]);
    }
  }

  DomainSegmentArray() = default;

  /// DomainSegmentArray destructor used to delete all pointers for
  /// clean-up
  ~DomainSegmentArray() {
    clear(); // deallocate all data
  }
};
} // namespace viennahrle

#endif // HRLE_DOMAIN_SEGMENT_ARRAY_HPP
