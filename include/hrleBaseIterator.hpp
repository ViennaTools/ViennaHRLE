#ifndef HRLE_BASE_ITERATOR_HPP
#define HRLE_BASE_ITERATOR_HPP

#include <iostream>

#include "hrleRunTypeValues.hpp"
#include "hrleSizeType.hpp"
#include "hrleIndexType.hpp"
#include "hrleVectorType.hpp"

/// Base class for RunsIterator and RunsIteratorOffset. Should/Can not be used
/// directly.
template <class hrleDomain> class hrleBaseIterator {
protected:
  static constexpr int D = hrleDomain::dimension;
  typedef typename std::conditional<std::is_const<hrleDomain>::value,
                                    const typename hrleDomain::hrleValueType,
                                    typename hrleDomain::hrleValueType>::type
      hrleValueType;

  hrleDomain &domain;
  hrleVectorType<hrleSizeType, D + 1> startIndicesPos;
  hrleVectorType<hrleSizeType, D> runTypePos;
  hrleVectorType<hrleIndexType, D> startRunAbsCoords;
  hrleVectorType<hrleIndexType, D> endRunAbsCoords;
  hrleVectorType<hrleIndexType, D> absCoords;
  hrleVectorType<hrleIndexType, D> end_absCoords;
  int r_level;
  int s_level;
  int sub;

  void go_up_BA() {
    ++r_level;
    // shfdhsfhdskjhgf assert(r_level==s_level);
  }

  void go_up_AB() {
    // shfdhsfhdskjhgf assert((r_level==s_level) || (s_level==r_level+1));
    s_level = r_level + 1;
    absCoords[r_level] = domain.getGrid().getMinIndex(r_level);
    end_absCoords[r_level] = domain.getGrid().getMaxIndex(r_level);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

public:
  void print() { // TODO
    std::cout << "startIndicesPos: " << startIndicesPos << std::endl;
    std::cout << "runTypePos: " << runTypePos << std::endl;
    std::cout << "startRunAbsCoords: " << startRunAbsCoords << std::endl;
    std::cout << "endRunAbsCoords: " << endRunAbsCoords << std::endl;
    std::cout << "absCoords: " << absCoords << std::endl;
    std::cout << "end_absCoords: " << end_absCoords << std::endl;
    std::cout << "r_level: " << r_level << std::endl;
    std::cout << "s_level: " << s_level << std::endl;
    std::cout << "sub: " << sub << std::endl;
    std::cout << "value: " << getValue() << std::endl;
    // std::cout << "sign: " << sign() << std::endl;
    // std::cout << "active: " << is_active() << std::endl;
    std::cout << "current segment: " << getSegmentId() << std::endl
              << std::endl;
  }

  hrleBaseIterator(hrleDomain &passedDomain)
      : domain(passedDomain), absCoords(domain.getGrid().getMinIndex()),
        end_absCoords(domain.getGrid().getMaxIndex()), r_level(D), s_level(D),
        sub(0) {
    startIndicesPos[D] = 0;
  }

  hrleDomain &getDomain() { return domain; }

  // const hrleBaseIterator &operator=(const hrleBaseIterator &it) {
  //   // copy assignment
  //   domain = it.domain;
  //   startIndicesPos = it.startIndicesPos;
  //   runTypePos = it.runTypePos;
  //   startRunAbsCoords = it.startRunAbsCoords;
  //   endRunAbsCoords = it.endRunAbsCoords;
  //   absCoords = it.absCoords;
  //   end_absCoords = it.end_absCoords;
  //   s_level = it.s_level;
  //   r_level = it.r_level;
  //   sub = it.sub;
  //   return *this;
  // }

  bool isFinished() const {
    // returns true if iterator reached the end
    return (r_level == D);
  }

  hrleIndexType getStartIndices(int dir) const {
    // returns the start index of a run for a certain axis direction
    return absCoords[dir];
  }

  const hrleVectorType<hrleIndexType, D> &getStartIndices() const {
    return absCoords;
  }

  hrleIndexType getEndIndices(int dir) const {
    // returns the end index of a run for a certain axis direction
    return end_absCoords[dir];
  }

  const hrleVectorType<hrleIndexType, D> &getEndIndices() const {
    return end_absCoords;
  }

  hrleSizeType getRunTypePosition() const {
    if (s_level == 0) {
      // shfdhsfhdskjhgf assert(s_level==r_level);
      return startIndicesPos[0];
    } else {
      // shfdhsfhdskjhgf assert(s_level>r_level);
      return runTypePos[r_level];
    }
  }

  int getLevel() const { return s_level; }

  hrleSizeType getRunCode() const { return startIndicesPos[r_level]; }

  hrleSizeType getSegmentId() const {
    // std::cout << "getSegmentId, r_level: " << r_level << std::endl;
    // std::cout << "getSegmentId, startIndicesPos[r_level]: " <<
    // startIndicesPos[r_level] << std::endl;
    const hrleSizeType r = getRunCode();
    // std::cout << "getSegmentId(): getRunCode: " << r << std::endl;
    // shfdhsfhdskjhgf assert(POS_PT<hrleRunTypeValues::SEGMENT_PT);
    // shfdhsfhdskjhgf assert(NEG_PT<hrleRunTypeValues::SEGMENT_PT);
    if (r >= hrleRunTypeValues::SEGMENT_PT) {
      // std::cout << "getSegmentId(): crc - SEGPT: " << r -
      // hrleRunTypeValues::SEGMENT_PT << std::endl;
      return r - hrleRunTypeValues::SEGMENT_PT;
    } else {
      return sub;
    }
  }

  // hrleSizeType getSegmentId() const { return sub; }

  hrleSizeType getPointId() const {
    // returns the getPointId if it is a defined run
    hrleSizeType tmp = getRunCode();
    if (isDefined())
      tmp += domain.pointIdOffsets[sub];
    return tmp;
  }

  bool isDefined() const {
    // returns if a run is defined or not
    // NOTE: if a run is defined, it always has length 1 and therefore it always
    //      represents a single grid point
    return (s_level == 0);
  }

  hrleValueType &getValue() {
    // returns the level set value for the current run
    // if the run is undefined either POS_VALUE or NEG_VALUE is returned
    if (isDefined()) {
      return getDefinedValue();
    } else {
      // assert(getRunCode()==NEG_PT); //TODO
      return domain.domainSegments[sub]
          .undefinedValues[hrleRunTypeValues::UNDEF_PT - getRunCode()];
    }
  }

  hrleValueType &getDefinedValue() {
    // the same as "value", however this function assumes
    // that the current run is defined, and therefore no check is required
    // if the run is undefined or not
    // shfdhsfhdskjhgf assert(getRunCode()>=0);
    // shfdhsfhdskjhgf
    // assert(getRunCode()<l.sub_levelsets[sub].distances.size());
    return domain.domainSegments[sub].definedValues[getRunCode()];
  }
};

// typedef for const iterator
template <class hrleDomain>
using hrleConstBaseIterator = hrleBaseIterator<const hrleDomain>;

#endif // HRLE_BASE_ITERATOR_HPP
