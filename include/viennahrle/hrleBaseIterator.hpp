#ifndef HRLE_BASE_ITERATOR_HPP
#define HRLE_BASE_ITERATOR_HPP

#include <cassert>
#include <iostream>

#include "hrleRunTypeValues.hpp"
#include "hrleTypes.hpp"

#include <vcVectorType.hpp>

/// Base class for RunsIterator and RunsIteratorOffset. Should/Can not be used
/// directly.
namespace viennahrle {
using namespace viennacore;
template <class hrleDomain> class BaseIterator {
protected:
  static constexpr int D = hrleDomain::dimension;
  typedef std::conditional_t<std::is_const_v<hrleDomain>,
                             const typename hrleDomain::ValueType,
                             typename hrleDomain::ValueType>
      ValueType;

  hrleDomain &domain;
  VectorType<SizeType, D + 1> startIndicesPos;
  VectorType<SizeType, D> runTypePos;
  VectorType<IndexType, D> startRunAbsCoords;
  VectorType<IndexType, D> endRunAbsCoords;
  VectorType<IndexType, D> absCoords;
  VectorType<IndexType, D> endAbsCoords;
  int r_level;
  int s_level;
  int sub;

  void go_up_BA() {
    ++r_level;
    assert(r_level == s_level);
  }

  void go_up_AB() {
    // assert((r_level == s_level) || (s_level == r_level + 1));
    s_level = r_level + 1;
    absCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);
    endAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);
  }

public:
  using DomainType = hrleDomain;

  void print() {
    std::cout << "startIndicesPos: " << startIndicesPos << std::endl;
    std::cout << "runTypePos: " << runTypePos << std::endl;
    std::cout << "startRunAbsCoords: " << startRunAbsCoords << std::endl;
    std::cout << "endRunAbsCoords: " << endRunAbsCoords << std::endl;
    std::cout << "absCoords: " << absCoords << std::endl;
    std::cout << "endAbsCoords: " << endAbsCoords << std::endl;
    std::cout << "r_level: " << r_level << std::endl;
    std::cout << "s_level: " << s_level << std::endl;
    std::cout << "sub: " << sub << std::endl;
    std::cout << "value: " << getValue() << std::endl;
    std::cout << "current segment: " << getSegmentRun() << std::endl
              << std::endl;
  }

  explicit BaseIterator(hrleDomain &passedDomain)
      : domain(passedDomain), absCoords(domain.getGrid().getMinGridPoint()),
        endAbsCoords(domain.getGrid().getMaxGridPoint()), r_level(D),
        s_level(D), sub(0) {
    startIndicesPos[D] = 0;
  }

  DomainType &getDomain() { return domain; }

  // const BaseIterator &operator=(const BaseIterator &it) {
  //   // copy assignment
  //   domain = it.domain;
  //   startIndicesPos = it.startIndicesPos;
  //   runTypePos = it.runTypePos;
  //   startRunAbsCoords = it.startRunAbsCoords;
  //   endRunAbsCoords = it.endRunAbsCoords;
  //   absCoords = it.absCoords;
  //   endAbsCoords = it.endAbsCoords;
  //   s_level = it.s_level;
  //   r_level = it.r_level;
  //   sub = it.sub;
  //   return *this;
  // }

  bool isFinished() const {
    // returns true if iterator reached the end
    return (r_level == D);
  }

  IndexType getStartIndices(int dir) const {
    // returns the start index of a run for a certain axis direction
    return absCoords[dir];
  }

  const VectorType<IndexType, D> &getStartIndices() const { return absCoords; }

  IndexType getEndIndices(int dir) const {
    // returns the end index of a run for a certain axis direction
    return endAbsCoords[dir];
  }

  const VectorType<IndexType, D> &getEndIndices() const { return endAbsCoords; }

  SizeType getRunTypePosition() const {
    if (s_level == 0) {
      assert(s_level == r_level);
      return startIndicesPos[0];
    } else {
      assert(s_level > r_level);
      return runTypePos[r_level];
    }
  }

  int getLevel() const { return s_level; }

  SizeType getRunCode() const { return startIndicesPos[r_level]; }

  SizeType getSegmentRun() const {
    const SizeType r = getRunCode();
    if (r >= RunTypeValues::SEGMENT_PT) {
      return r - RunTypeValues::SEGMENT_PT;
    } else {
      return sub;
    }
  }

  SizeType getSegmentId() const { return sub; }

  SizeType getPointId() const {
    // returns the getPointId if it is a defined run
    SizeType tmp = getRunCode();
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

  const ValueType &getValue() const {
    // returns the level set value for the current run
    // if the run is undefined either POS_VALUE or NEG_VALUE is returned
    if (isDefined()) {
      return getDefinedValue();
    } else {
      // assert(getRunCode()==NEG_PT); //TODO
      return domain.domainSegments[sub]
          .undefinedValues[getRunCode() - RunTypeValues::UNDEF_PT];
    }
  }

  ValueType &getValue() {
    return const_cast<ValueType &>(
        const_cast<const BaseIterator *>(this)->getValue());
  }

  const ValueType &getDefinedValue() const {
    // the same as "value", however this function assumes
    // that the current run is defined, and therefore no check is required
    // if the run is undefined or not
    // assert(getRunCode() < l.sub_levelsets[sub].distances.size());
    return domain.domainSegments[sub].definedValues[getRunCode()];
  }

  ValueType &getDefinedValue() {
    return const_cast<ValueType &>(
        const_cast<const BaseIterator *>(this)->getDefinedValue());
  }
};

// typedef for const iterator
template <class hrleDomain>
using ConstBaseIterator = BaseIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_BASE_ITERATOR_HPP
