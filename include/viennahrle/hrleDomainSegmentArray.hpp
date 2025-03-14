#ifndef HRLE_DOMAIN_SEGMENT_ARRAY_HPP
#define HRLE_DOMAIN_SEGMENT_ARRAY_HPP

#include "hrleDomainSegment.hpp"

template <class T = double, int D = 3> class hrleDomainSegmentArray {
  // contains the vector of hrleDomainSegment items that are used for
  // parallelization

  typedef hrleDomainSegment<T, D>
      *segmentPointer; // pointer to memory containing
                       // hrleDomainSegment object

  typedef std::vector<segmentPointer> hrleDomainSegmentPointerArray;
  hrleDomainSegmentPointerArray domainSegmentPointers;

public:
  // assignment operator and copy-constructor are deleted since vector only
  // holds pointers, use deepCopy or shallowCopy instead
  hrleDomainSegmentArray &operator=(const hrleDomainSegmentArray &s) = delete;
  hrleDomainSegmentArray(const hrleDomainSegmentArray &s) = delete;

  typedef typename hrleDomainSegmentPointerArray::size_type size_type;

  const hrleDomainSegment<T, D> &operator[](size_type i) const {
    return *(domainSegmentPointers[i]);
  }

  hrleDomainSegment<T, D> &operator[](size_type i) {
    return *(domainSegmentPointers[i]);
  }

  const hrleDomainSegment<T, D> &back() const {
    return *(domainSegmentPointers.back());
  }

  hrleDomainSegment<T, D> &back() { return *(domainSegmentPointers.back()); }

  void swap(hrleDomainSegmentArray &s) noexcept {
    domainSegmentPointers.swap(s.domainSegmentPointers);
  }

  void clear() {
    for (typename hrleDomainSegmentPointerArray::iterator it =
             domainSegmentPointers.begin();
         it != domainSegmentPointers.end(); ++it)
      if (*it != 0)
        delete *it;
    domainSegmentPointers.clear();
  }

  /// initialize the parallelization by implementing i number of
  /// hrleDomainSegment objects through a segmentPointer object, that
  /// each correspond to a parallel core
  void initialize(size_type i, const hrleGrid<D> &g,
                  hrleAllocationType<hrleSizeType, D> a =
                      hrleAllocationType<hrleSizeType, D>()) {
    clear(); // delete all data in the array

    domainSegmentPointers.resize(i);

    a *= a.allocationFactor; // increase allocation by extra margin defined in
                             // hrleAllocationType
    a /= i;

#pragma omp parallel for schedule(static, 1)
    // parallelization - Iterations divided into chunks of size 1.
    // Each chunk is assigned to a thread
    for (int k = 0; k < static_cast<int>(domainSegmentPointers.size()); ++k) {
      domainSegmentPointers[k] =
          segmentPointer(new hrleDomainSegment<T, D>(g, a));
    }
  }

  size_type getNumberOfSegments() const {
    // return the size of the vector for parallelization
    return domainSegmentPointers.size();
  }

  void deepCopy(const hrleGrid<D> *g, const hrleDomainSegmentArray &s) {
    clear(); // deallocate all domain segments
    domainSegmentPointers.resize(s.getNumberOfSegments());

#pragma omp parallel for schedule(static, 1)
    // parallelization - Iterations divided into chunks of size 1.
    // Each chunk is assigned to a thread
    for (int k = 0; k < static_cast<int>(s.getNumberOfSegments()); ++k) {
      domainSegmentPointers[k] =
          segmentPointer(new hrleDomainSegment<T, D>(*g, s[k]));
    }
  }

  void shallowCopy(const hrleDomainSegmentArray &s) {
    for (int k = 0; k < static_cast<int>(s.getNumberOfSegments()); ++k) {
      domainSegmentPointers[k] = segmentPointer(&s[k]);
    }
  }

  hrleDomainSegmentArray() = default;

  /// hrleDomainSegmentArray destructor used to delete all pointers for
  /// clean-up
  ~hrleDomainSegmentArray() {
    clear(); // deallocate all data
  }
};

#endif // HRLE_DOMAIN_SEGMENT_ARRAY_HPP
