#ifndef HRLE_RUNS_ITERATOR_HPP
#define HRLE_RUNS_ITERATOR_HPP

#include "hrleBaseIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// This iterator iterates over all runs and stops at every start index of each
/// run. A run can be defined or undefined: A defined run
/// corresponds to a single defined grid point.
template <class hrleDomain>
class SparseIterator : public BaseIterator<hrleDomain> {

  typedef typename hrleDomain::DomainSegmentType DomainSegmentType;

  using BaseIterator<hrleDomain>::domain;
  using BaseIterator<hrleDomain>::r_level;
  using BaseIterator<hrleDomain>::s_level;
  using BaseIterator<hrleDomain>::sub;
  using BaseIterator<hrleDomain>::absCoords;
  using BaseIterator<hrleDomain>::endAbsCoords;
  using BaseIterator<hrleDomain>::startRunAbsCoords;
  using BaseIterator<hrleDomain>::endRunAbsCoords;
  using BaseIterator<hrleDomain>::runTypePos;
  using BaseIterator<hrleDomain>::startIndicesPos;
  using BaseIterator<hrleDomain>::go_up_AB;
  using BaseIterator<hrleDomain>::go_up_BA;
  using BaseIterator<hrleDomain>::D;

private:
  // Cache domain segment reference for performance
  mutable const DomainSegmentType *currentSegment;

  // Update cached segment when sub changes
  void updateSegmentCache() const noexcept {
    currentSegment = &domain.domainSegments[sub];
  }

  // Get cached domain segment
  const DomainSegmentType &getDomainSegment() const noexcept {
    return *currentSegment;
  }

  void go_down_AB(IndexType c) { // find right runTypePos
    const DomainSegmentType &sl = getDomainSegment();
    const SizeType s = startIndicesPos[s_level]; // Remove unnecessary reference

    --r_level;                                // go down one level
    SizeType r = sl.startIndices[r_level][s]; // r is the start index of the run

    auto start_breaks = sl.runBreaks[r_level].begin() + (r - s);
    auto end_breaks = sl.runBreaks[r_level].begin() +
                      (sl.getStartIndex(r_level, s + 1) - (s + 1));

    assert(start_breaks <= end_breaks);

    auto pos_breaks = std::upper_bound(start_breaks, end_breaks, c);

    r += static_cast<SizeType>(pos_breaks - start_breaks);

    if (pos_breaks == start_breaks) {
      // position c is before the first break -> position is undefined
      startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);
    } else {
      // position c is between two breaks -> position is in defined run
      assert(pos_breaks > start_breaks);
      startRunAbsCoords[r_level] = *(pos_breaks - 1);
    }

    if (pos_breaks == end_breaks) {
      // position c is after the last break -> position is undefined
      endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);
    } else {
      // position c is between two breaks -> position is in defined run
      // start run = end run -> single defined run
      endRunAbsCoords[r_level] = (*pos_breaks) - 1;
    }

    runTypePos[r_level] = r;

    assert(s_level >= 1);
    assert(s_level == r_level + 1);
  }

  void go_down_AB_first() { // find right runTypePos
    assert(s_level == r_level);

    const DomainSegmentType &sl = getDomainSegment();
    const SizeType s = startIndicesPos[s_level]; // Remove unnecessary reference

    --r_level;

    const SizeType r = sl.startIndices[r_level][s];

    startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);

    endRunAbsCoords[r_level] = sl.getRunEndCoord(r_level, s, r);

    runTypePos[r_level] = r;

    assert(s_level >= 1);
    assert(s_level == r_level + 1);
  }

  void go_down_AB_last() { // find right runTypePos
    assert(s_level == r_level);

    const DomainSegmentType &sl = getDomainSegment();
    const SizeType s = startIndicesPos[s_level]; // Remove unnecessary reference

    --r_level;

    const SizeType r = sl.getStartIndex(r_level, s + 1) - 1;

    startRunAbsCoords[r_level] = sl.getRunStartCoord(r_level, s, r);

    endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);

    runTypePos[r_level] = r;

    assert(s_level >= 1);
    assert(s_level == r_level + 1);
  }

  bool go_down_BA(IndexType c) {
    assert(s_level == r_level + 1);
    assert(s_level >= 1);

    const DomainSegmentType &sl = getDomainSegment();

    startIndicesPos[s_level - 1] = sl.runTypes[r_level][runTypePos[r_level]];
    if (DomainSegmentType::isPtIdDefined(startIndicesPos[s_level - 1])) {
      --s_level;
      startIndicesPos[s_level] += (c - startRunAbsCoords[r_level]);
      absCoords[s_level] = c;
      endAbsCoords[s_level] = c;
      return true;
    } else {
      absCoords[r_level] = startRunAbsCoords[r_level];
      endAbsCoords[r_level] = endRunAbsCoords[r_level];
      return false;
    }
  }

  bool go_next_A() noexcept {
    if (s_level == r_level) {
      // if
      // (startIndicesPos[s_level]<l.runtypes[r_level][runTypePos[r_level]]+(end_run_coords[r_level]-start_run_coords[r_level]))
      // {
      if (absCoords[s_level] < endRunAbsCoords[r_level]) {
        ++startIndicesPos[s_level];
        ++absCoords[s_level];
        ++endAbsCoords[s_level];
        return true;
      }
    }
    return false;
  }

  bool go_next_B() {
    assert(s_level == r_level + 1);
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxIndex(r_level));

    const DomainSegmentType &sl = getDomainSegment();

    if (endRunAbsCoords[r_level] != domain.getGrid().getMaxGridPoint(r_level)) {
      // current run is not the last run in the segment
      ++runTypePos[r_level];
      startRunAbsCoords[r_level] = endRunAbsCoords[r_level] + 1;
      endRunAbsCoords[r_level] = sl.getRunEndCoord(
          r_level, startIndicesPos[s_level], runTypePos[r_level]);
      return true;
    } else {
      return false;
    }
  }

  bool go_previous_A() noexcept {
    if (s_level == r_level) {
      if (absCoords[s_level] > startRunAbsCoords[r_level]) {
        --startIndicesPos[s_level];
        --absCoords[s_level];
        --endAbsCoords[s_level];
        return true;
      }
    }
    return false;
  }

  bool go_previous_B() {
    assert(s_level == r_level + 1);
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinIndex(r_level));

    const DomainSegmentType &sl = getDomainSegment();

    if (startRunAbsCoords[r_level] !=
        domain.getGrid().getMinGridPoint(r_level)) {
      --runTypePos[r_level];
      endRunAbsCoords[r_level] = startRunAbsCoords[r_level] - 1;
      startRunAbsCoords[r_level] = sl.getRunStartCoord(
          r_level, startIndicesPos[s_level], runTypePos[r_level]);
      return true;
    } else {
      return false;
    }
  }

public:
  using DomainType = hrleDomain;

  void print() { BaseIterator<hrleDomain>::print(); }

  explicit SparseIterator(hrleDomain &passedDomain, bool reverse = false)
      : BaseIterator<hrleDomain>(passedDomain), currentSegment(nullptr) {
    updateSegmentCache(); // Initialize cache
    if (reverse) {
      goToIndices(domain.getGrid().getMaxGridPoint());
    } else {
      goToIndices(domain.getGrid().getMinGridPoint());
    }
  }

  template <class V>
  SparseIterator(hrleDomain &lx, const V &v)
      : BaseIterator<hrleDomain>(lx), currentSegment(nullptr) {
    updateSegmentCache(); // Initialize cache
    goToIndices(v);
  }

  // TODO: it would make more sense if the call to next returned a bool whether
  // it was possible to move or not. currently calling next on a finished
  // iterator results in a segfault because r_level==D in
  // endRunAbsCoords[r_level] in go_next_B, it would make more sense to check
  // finished first, do nothing if it is and return false
  SparseIterator &operator++() {
    while (true) {
      if (go_next_A())
        break;
      go_up_AB();
      if (go_next_B()) {
        go_down_BA(startRunAbsCoords[r_level]);
        break;
      }
      go_up_BA();
      // shfdhsfhdskjhgf assert(it.r_level==it.s_level);
      if (r_level == D) {
        absCoords = domain.getGrid().incrementIndices(
            domain.getGrid().getMaxGridPoint());
        endAbsCoords = absCoords; // TODO
        // it.startIndicesPos[D]=it.l.num_pts();
        return *this;
      }
    }

    while (r_level == s_level && s_level > 0) {
      const IndexType c = domain.getGrid().getMinGridPoint(r_level - 1);
      // go_down_AB(c);
      go_down_AB_first();
      go_down_BA(c);
    }

    // check if new point is in the same segment, if not change segments
    const int newSub = BaseIterator<hrleDomain>::getSegmentRun();
    if (newSub != sub) {
      sub = newSub;
      updateSegmentCache(); // Update cache when segment changes
      goToIndices(newSub, absCoords);
    }

    return *this;
  }

  SparseIterator operator++(int) {
    SparseIterator temp(*this);
    operator++();
    return temp;
  }

  SparseIterator &operator--() {
    while (true) {
      if (go_previous_A())
        break;
      go_up_AB();
      if (go_previous_B()) {
        go_down_BA(endRunAbsCoords[r_level]);
        break;
      }
      go_up_BA();
      // shfdhsfhdskjhgf assert(it.r_level==it.s_level);
      if (r_level == D) {
        absCoords = domain.getGrid().decrementIndices(
            domain.getGrid().getMinGridPoint());
        endAbsCoords = absCoords; // TODO
        return *this;
      }
    }

    while (r_level == s_level && s_level > 0) {
      const IndexType c = domain.getGrid().getMaxGridPoint(r_level - 1);
      // go_down_AB(c);
      go_down_AB_last();
      go_down_BA(c);
    }

    const int newSub = BaseIterator<hrleDomain>::getSegmentRun();
    if (newSub != sub) {
      sub = newSub;
      updateSegmentCache(); // Update cache when segment changes
      goToIndices(newSub, BaseIterator<hrleDomain>::getEndIndices());
    }
    return *this;
  }

  SparseIterator operator--(int) {
    SparseIterator temp(*this);
    operator--();
    return temp;
  }

  /// safe version of operator++(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool next() {
    if (this->isFinished())
      return false;
    operator++();
    return true;
  }

  /// safe version of operator--(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool previous() {
    if (this->isFinished())
      return false;
    operator--();
    return true;
  }

  /// Advances the iterator to position v.
  /// If v is lexicographically higher than the current position
  /// the iterator will be moved back to v.
  /// If v is lexicographically smaller than the current position
  /// then the iterator will be moved until it reaches v
  template <class V> void goToIndicesSequential(const V &v) {
    if (v >= absCoords) {
      while (v > endAbsCoords) {
        ++(*this);
      }
    } else {
      while (v < absCoords) {
        --(*this);
      }
    }
  }

  template <class V> void goToIndices(const V &v) {
    goToIndices(0, v); // TODO
    const int newSub = this->getSegmentRun();
    if (newSub != 0) {
      sub = newSub;
      updateSegmentCache(); // Update cache when segment changes
      goToIndices(newSub, v);
    }
  }

  template <class V> void goToIndices(int subDomain, const V &v) {
    r_level = D;
    s_level = D;
    if (sub != subDomain) {
      sub = subDomain;
      updateSegmentCache(); // Update cache when segment changes
    }
    startIndicesPos[D] = 0;
    do {
      const IndexType c = v[r_level - 1];

      // shfdhsfhdskjhgf assert(isDefined(it.startIndicesPos[it.s_level]));

      go_down_AB(c);

      // shfdhsfhdskjhgf assert(c>=it.startRunAbsCoords[it.r_level]);
      // shfdhsfhdskjhgf assert(it.endRunAbsCoords[it.r_level]>=c);

      go_down_BA(c);
    } while (r_level == s_level && s_level > 0);
    // shfdhsfhdskjhgf assert(!it.isFinished());

    for (int h = 0; h < r_level; ++h) {
      absCoords[h] = domain.getGrid().getMinGridPoint(h);
      endAbsCoords[h] = domain.getGrid().getMaxGridPoint(h);
    }
  }
};

// typedef for const iterator
template <class hrleDomain>
using ConstSparseIterator = SparseIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_RUNS_ITERATOR_HPP
