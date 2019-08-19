#ifndef HRLE_RUNS_ITERATOR_HPP
#define HRLE_RUNS_ITERATOR_HPP

#include "hrleBaseIterator.hpp"
#include "hrleIndexType.hpp"

/// This iterator iterates over all runs and stops at every start index of each
/// run. A run can be defined or undefined: A defined run
/// corresponds to a single defined grid point.
template <class hrleDomain>
class hrleRunsIterator : public hrleBaseIterator<hrleDomain> {

  typedef typename hrleDomain::hrleDomainSegmentType hrleDomainSegmentType;

  using hrleBaseIterator<hrleDomain>::domain;
  using hrleBaseIterator<hrleDomain>::r_level;
  using hrleBaseIterator<hrleDomain>::s_level;
  using hrleBaseIterator<hrleDomain>::sub;
  using hrleBaseIterator<hrleDomain>::absCoords;
  using hrleBaseIterator<hrleDomain>::endAbsCoords;
  using hrleBaseIterator<hrleDomain>::startRunAbsCoords;
  using hrleBaseIterator<hrleDomain>::endRunAbsCoords;
  using hrleBaseIterator<hrleDomain>::runTypePos;
  using hrleBaseIterator<hrleDomain>::startIndicesPos;
  using hrleBaseIterator<hrleDomain>::go_up_AB;
  using hrleBaseIterator<hrleDomain>::go_up_BA;
  using hrleBaseIterator<hrleDomain>::D;

private:
  void go_down_AB(hrleIndexType c) { // find right runTypePos
    const hrleDomainSegmentType &sl = domain.domainSegments[sub];
    const hrleSizeType &s = startIndicesPos[s_level];

    --r_level;
    hrleSizeType r = sl.startIndices[r_level][s];

    typename std::vector<hrleIndexType>::const_iterator start_breaks =
        sl.runBreaks[r_level].begin() + (r - s);
    typename std::vector<hrleIndexType>::const_iterator end_breaks =
        sl.runBreaks[r_level].begin() +
        (sl.getStartIndex(r_level, s + 1) - (s + 1));

    // shfdhsfhdskjhgf assert(start_breaks<=end_breaks);

    typename std::vector<hrleIndexType>::const_iterator pos_breaks =
        std::upper_bound(start_breaks, end_breaks, c);

    r += hrleSizeType(pos_breaks - start_breaks);

    if (pos_breaks == start_breaks) {
      startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);
    } else {
      // shfdhsfhdskjhgf assert(pos_breaks>start_breaks);
      startRunAbsCoords[r_level] = *(pos_breaks - 1);
    }

    if (pos_breaks == end_breaks) {
      endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);
    } else {
      endRunAbsCoords[r_level] = (*pos_breaks) - 1;
    }

    runTypePos[r_level] = r;

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

  void go_down_AB_first() { // find right runTypePos

    // shfdhsfhdskjhgf assert(s_level==r_level);

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    const hrleSizeType &s = startIndicesPos[s_level];

    --r_level;

    const hrleSizeType r = sl.startIndices[r_level][s];

    startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);

    endRunAbsCoords[r_level] = sl.getRunEndCoord(r_level, s, r);

    runTypePos[r_level] = r;

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

  void go_down_AB_last() { // find right runTypePos

    // shfdhsfhdskjhgf assert(s_level==r_level);

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    const hrleSizeType &s = startIndicesPos[s_level];

    --r_level;

    const hrleSizeType r = sl.getStartIndex(r_level, s + 1) - 1;

    startRunAbsCoords[r_level] = sl.getRunStartCoord(r_level, s, r);

    endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);

    runTypePos[r_level] = r;

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

  bool go_down_BA(hrleIndexType c) {
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf assert(s_level>=1);
    // std::cout << "go_down_BA:" << std::endl;
    // print();

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    startIndicesPos[s_level - 1] = sl.runTypes[r_level][runTypePos[r_level]];
    if (hrleDomainSegmentType::isPtIdDefined(startIndicesPos[s_level - 1])) {
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

  bool go_next_A() {
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
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxIndex(r_level));

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    if (endRunAbsCoords[r_level] != domain.getGrid().getMaxGridPoint(r_level)) {
      ++runTypePos[r_level];
      startRunAbsCoords[r_level] = endRunAbsCoords[r_level] + 1;
      endRunAbsCoords[r_level] = sl.getRunEndCoord(
          r_level, startIndicesPos[s_level], runTypePos[r_level]);
      return true;
    } else {
      return false;
    }
  }

  bool go_previous_A() {
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
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinIndex(r_level));

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

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
  void print() { hrleBaseIterator<hrleDomain>::print(); }

  hrleRunsIterator(hrleDomain &passedDomain, bool reverse = false)
      : hrleBaseIterator<hrleDomain>(passedDomain) {
    if (reverse) {
      goToIndices(domain.getGrid().getMaxGridPoint());
    } else {
      goToIndices(domain.getGrid().getMinGridPoint());
    }
  }

  template <class V>
  hrleRunsIterator(hrleDomain &lx, const V &v)
      : hrleBaseIterator<hrleDomain>(lx) {
    goToIndices(v);
  }

  // TODO: it would make more sense if the call to next returned a bool whether
  // it was possible to move or not. currently calling next on a finished
  // iterator results in a segfault because r_level==D in
  // endRunAbsCoords[r_level] in go_next_B, it would make more sense to check
  // finished first, do nothing if it is and return false
  hrleRunsIterator<hrleDomain> &operator++() {
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
        absCoords =
            domain.getGrid().incrementIndices(domain.getGrid().getMaxIndex());
        endAbsCoords = absCoords; // TODO
        // it.startIndicesPos[D]=it.l.num_pts();
        return *this;
      }
    }

    while (r_level == s_level && s_level > 0) {
      const hrleIndexType c = domain.getGrid().getMinIndex(r_level - 1);
      // go_down_AB(c);
      go_down_AB_first();
      go_down_BA(c);
    }

    // check if new point is in the same segment, if not change segments
    int s = hrleBaseIterator<hrleDomain>::getSegmentId();
    if (s != sub) {
      goToIndices(s, absCoords);
    }

    return *this;
  }

  hrleRunsIterator<hrleDomain> operator++(int) {
    hrleRunsIterator temp(*this);
    operator++();
    return temp;
  }

  hrleRunsIterator<hrleDomain> &operator--() {
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
        absCoords =
            domain.getGrid().decrementIndices(domain.getGrid().getMinIndex());
        endAbsCoords = absCoords; // TODO
        return *this;
      }
    }

    while (r_level == s_level && s_level > 0) {
      const hrleIndexType c = domain.getGrid().getMaxIndex(r_level - 1);
      // go_down_AB(c);
      go_down_AB_last();
      go_down_BA(c);
    }

    int s = hrleBaseIterator<hrleDomain>::getSegmentId();
    if (s != sub) {
      goToIndices(s, hrleBaseIterator<hrleDomain>::getEndIndices());
      /*if (s<it.l.segmentation.size()) {
          //shfdhsfhdskjhgf
      assert(it.l.grid().decrementIndices(it.l.segmentation[s])==it.getEndIndices());
      } else {
          //shfdhsfhdskjhgf
      assert(it.l.grid().max_point_index()==it.getEndIndices());
      }*/
    }
    return *this;
  }

  hrleRunsIterator<hrleDomain> operator--(int) {
    hrleRunsIterator temp(*this);
    operator--();
    return temp;
  }

  /// safe version of operator++(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool next() {
    if (hrleBaseIterator<hrleDomain>::isFinished())
      return false;
    operator++();
    return true;
  }

  /// safe version of operator--(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool previous() {
    if (hrleBaseIterator<hrleDomain>::isFinished())
      return false;
    operator--();
    return true;
  }

  template <class V> void goToIndices(const V &v) {
    goToIndices(0, v); // TODO
    int s = hrleBaseIterator<hrleDomain>::getSegmentId();
    if (s != 0)
      goToIndices(s, v);
  }

  template <class V> void goToIndices(int subDomain, const V &v) {
    r_level = D;
    s_level = D;
    sub = subDomain;
    startIndicesPos[D] = 0;
    do {
      const hrleIndexType c = v[r_level - 1];

      // shfdhsfhdskjhgf assert(isDefined(it.startIndicesPos[it.s_level]));

      go_down_AB(c);

      // shfdhsfhdskjhgf assert(c>=it.startRunAbsCoords[it.r_level]);
      // shfdhsfhdskjhgf assert(it.endRunAbsCoords[it.r_level]>=c);

      go_down_BA(c);
    } while (r_level == s_level && s_level > 0);
    // shfdhsfhdskjhgf assert(!it.isFinished());

    for (int h = 0; h < r_level; ++h) {
      absCoords[h] = domain.getGrid().getMinIndex(h);
      endAbsCoords[h] = domain.getGrid().getMaxIndex(h);
    }
  }
};

// typedef for const iterator
template <class hrleDomain>
using hrleConstRunsIterator = hrleRunsIterator<const hrleDomain>;

#endif // HRLE_RUNS_ITERATOR_HPP
