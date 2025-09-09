#ifndef HRLE_OFFSET_RUNS_ITERATOR_HPP
#define HRLE_OFFSET_RUNS_ITERATOR_HPP

#include <bitset>
#include <cassert>
#include <vector>

#include "hrleBaseIterator.hpp"

namespace viennahrle {
using namespace viennacore;
/// The SparseOffsetIterator iterates over the runs exactly like
/// hrleSparseIterator, but is set at indices given by the offset it was
/// initialised with. So if the offset was (0,0,1) and this iterator is
/// currently at the point (2,2,2), calling getValue() will return the value at
/// (2,2,2), but calling indices() will return (2,2,3). By comparing its indices
/// with those of an hrleSparseIterator and only incrementing when necessary,
/// this iterator can be used to keep track of the cartesian neighbour of an
/// hrleSparseIterator, taking into account boundary conditions.
template <class hrleDomain>
class SparseOffsetIterator : public BaseIterator<hrleDomain> {

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

  Index<D> start_run_rel_coords;
  Index<D> end_run_rel_coords;
  Index<D> rel_coords;
  Index<D> offset;
  std::bitset<D> move_inverse;

protected:
  void go_down_AB(IndexType abs_c) { // X         //find right runTypePos

    // assert(s_level==r_level);
    // assert(domain.isDefined(startIndicesPos[s_level]));

    const DomainSegmentType &sl = domain.domainSegments[sub];

    const auto &s = startIndicesPos[s_level];

    --r_level;

    int cycles = 0;
    IndexType rel_c = domain.getGrid().globalIndex2LocalIndex(
        r_level, abs_c, offset[r_level], cycles);

    // assert(s<sl.startIndices[r_level].size());
    // assert(r_level>=0);

    auto r = sl.startIndices[r_level][s];

    auto start_breaks = sl.runBreaks[r_level].begin() + (r - s);
    auto end_breaks = sl.runBreaks[r_level].begin() +
                      (sl.getStartIndex(r_level, s + 1) - (s + 1));

    assert(start_breaks <= end_breaks);

    auto pos_breaks = std::upper_bound(start_breaks, end_breaks, rel_c);

    r += static_cast<SizeType>(pos_breaks - start_breaks);

    runTypePos[r_level] = r;

    if (pos_breaks == start_breaks) {
      start_run_rel_coords[r_level] = domain.getGrid().getMinGridPoint(r_level);
    } else {
      // assert(pos_breaks>start_breaks);
      start_run_rel_coords[r_level] = *(pos_breaks - 1);
    }

    if (pos_breaks == end_breaks) {
      end_run_rel_coords[r_level] = domain.getGrid().getMaxGridPoint(r_level);
    } else {
      // assert(pos_breaks<end_breaks);
      end_run_rel_coords[r_level] = (*pos_breaks) - 1;
    }

    // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
    // assert(start_run_rel_coords[r_level]<=rel_c);
    // assert(rel_c<=end_run_rel_coords[r_level]);
    // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

    // calculate run_absCoords
    if (domain.getGrid().isBoundaryPeriodic(r_level)) {

      const IndexType &rel_s = start_run_rel_coords[r_level];
      const IndexType &rel_e = end_run_rel_coords[r_level];

      move_inverse.reset(r_level);

      startRunAbsCoords[r_level] =
          std::max(domain.getGrid().localIndex2GlobalIndex(
                       r_level, rel_s, cycles, offset[r_level]),
                   domain.getGrid().getMinGridPoint(r_level));
      endRunAbsCoords[r_level] =
          std::min(domain.getGrid().localIndex2GlobalIndex(
                       r_level, rel_e, cycles, offset[r_level]),
                   domain.getGrid().getMaxGridPoint(r_level));

    } else {

      move_inverse.set(r_level, cycles & 1);

      if (start_breaks == end_breaks) {

        startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);
        endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);

      } else {

        IndexType rel_s = start_run_rel_coords[r_level];
        IndexType rel_e = end_run_rel_coords[r_level];

        if (pos_breaks == start_breaks)
          rel_s +=
              (start_run_rel_coords[r_level] - end_run_rel_coords[r_level]);
        else if (pos_breaks == end_breaks)
          rel_e +=
              (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);

        if ((cycles & 1) != 0)
          std::swap(rel_s, rel_e);

        startRunAbsCoords[r_level] =
            std::max(domain.getGrid().localIndex2GlobalIndex(
                         r_level, rel_s, cycles, offset[r_level]),
                     domain.getGrid().getMinGridPoint(r_level));

        endRunAbsCoords[r_level] =
            std::min(domain.getGrid().localIndex2GlobalIndex(
                         r_level, rel_e, cycles, offset[r_level]),
                     domain.getGrid().getMaxGridPoint(r_level));
      }
    }

    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
    // assert(startRunAbsCoords[r_level]<=abs_c);
    // assert(abs_c<=endRunAbsCoords[r_level]);
    // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);
    // assert(s_level>=1);
    // assert(s_level==r_level+1);

    if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
        end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
      // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
      // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
      // ||
      // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
    }

    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
    // offset[r_level])>=start_run_rel_coords[r_level]);
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
    // offset[r_level])<=end_run_rel_coords[r_level]);

    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
    // offset[r_level])>=start_run_rel_coords[r_level]);
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
    // offset[r_level])<=end_run_rel_coords[r_level]);
  }

  void go_down_AB_first() { // find right runTypePos

    // assert(s_level==r_level);

    const IndexType c = domain.getGrid().getMinGridPoint(r_level - 1);
    // assert(l.isDefined(startIndicesPos[s_level]));
    go_down_AB(c);

    // assert(s_level>=1);
    // assert(s_level==r_level+1);
  }

  void go_down_AB_last() { // find right runTypePos

    // assert(s_level==r_level);

    const IndexType c = domain.getGrid().getMaxGridPoint(r_level - 1);
    // assert(l.isDefined(startIndicesPos[s_level]));
    go_down_AB(c);

    // assert(s_level>=1);
    // assert(s_level==r_level+1);
  }

  bool go_down_BA(IndexType abs_c) { // X

    // assert(s_level==r_level+1);
    // assert(s_level>=1);

    const DomainSegmentType &sl = domain.domainSegments[sub];

    startIndicesPos[s_level - 1] = sl.runTypes[r_level][runTypePos[r_level]];
    if (sl.isPtIdDefined(startIndicesPos[s_level - 1])) {
      --s_level;
      // assert(domain.isDefined(startIndicesPos[s_level]));

      IndexType rel_c = domain.getGrid().globalIndex2LocalIndex(
          r_level, abs_c, offset[r_level]);

      // assert(rel_c>=start_run_rel_coords[r_level]);
      // assert(rel_c<=end_run_rel_coords[r_level]);

      startIndicesPos[s_level] += (rel_c - start_run_rel_coords[r_level]);
      rel_coords[s_level] = rel_c;
      absCoords[s_level] = abs_c;
      endAbsCoords[s_level] = abs_c;

      // assert(domain.isDefined(startIndicesPos[s_level]));

      return true;
    } else {
      absCoords[r_level] = startRunAbsCoords[r_level];
      endAbsCoords[r_level] = endRunAbsCoords[r_level];

      // assert(s_level==r_level+1);
      return false;
    }
  }

  bool go_next_A() { // X

    if (s_level == r_level) {

      // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
      // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
      // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
      // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

      // if
      // (startIndicesPos[s_level]<domain.runtypes[r_level][runTypePos[r_level]]+(end_run_coords[r_level]-start_run_coords[r_level]))
      // {

      if (absCoords[s_level] < endRunAbsCoords[r_level]) {
        ++absCoords[s_level];
        ++endAbsCoords[s_level];

        if (domain.getGrid().isBoundaryPeriodic(r_level)) {
          // assert(rel_coords[s_level]!=domain.getGrid().getMaxGridPoint(s_level));
          ++rel_coords[s_level];
          ++startIndicesPos[s_level];
        } else {

          if (rel_coords[s_level] ==
              domain.getGrid().getMaxGridPoint(s_level)) {
            move_inverse.set(s_level);
          } else if (rel_coords[s_level] ==
                     domain.getGrid().getMinGridPoint(s_level)) {
            move_inverse.reset(s_level);
          }

          if (move_inverse.test(s_level)) {
            // assert(startIndicesPos[s_level]>0);
            --startIndicesPos[s_level];
            --rel_coords[s_level];
          } else {
            ++startIndicesPos[s_level];
            ++rel_coords[s_level];
          }
        }

        // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
        // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
        // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
        // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

        return true;
      }
    }
    return false;
  }

  bool go_previous_A() {

    if (s_level == r_level) {

      // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
      // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
      // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
      // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

      if (absCoords[s_level] > startRunAbsCoords[r_level]) {
        --absCoords[s_level];
        --endAbsCoords[s_level];

        if (domain.getGrid().isBoundaryPeriodic(r_level)) {
          // assert(rel_coords[s_level]!=domain.getGrid().getMinGridPoint(s_level));
          --rel_coords[s_level];
          --startIndicesPos[s_level];
        } else {
          if (rel_coords[s_level] ==
              domain.getGrid().getMaxGridPoint(s_level)) {
            move_inverse.reset(s_level);
          } else if (rel_coords[s_level] ==
                     domain.getGrid().getMinGridPoint(s_level)) {
            move_inverse.set(s_level);
          }
          if (move_inverse.test(s_level)) {
            ++startIndicesPos[s_level];
            ++rel_coords[s_level];
          } else {
            --startIndicesPos[s_level];
            --rel_coords[s_level];
          }
        }

        // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
        // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
        // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
        // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

        return true;
      }
    }
    return false;
  }

  bool go_next_B() {

    // assert(s_level==r_level+1);
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxGridPoint(r_level));
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));

    const DomainSegmentType &sl = domain.domainSegments[sub];

    // assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
    // startIndicesPos[s_level]));
    // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
    // startIndicesPos[s_level]+1));

    if (endRunAbsCoords[r_level] != domain.getGrid().getMaxGridPoint(r_level)) {

      if (domain.getGrid().isBoundaryPeriodic(
              r_level)) { // periodic boundary conditions

        if (end_run_rel_coords[r_level] ==
            domain.getGrid().getMaxGridPoint(r_level)) {
          runTypePos[r_level] =
              sl.getStartIndex(r_level, startIndicesPos[s_level]);
          start_run_rel_coords[r_level] =
              domain.getGrid().getMinGridPoint(r_level);
        } else {
          ++runTypePos[r_level];
          start_run_rel_coords[r_level] = end_run_rel_coords[r_level] + 1;
        }

        end_run_rel_coords[r_level] = sl.getRunEndCoord(
            r_level, startIndicesPos[s_level], runTypePos[r_level]);

        ++endRunAbsCoords[r_level];
        startRunAbsCoords[r_level] = endRunAbsCoords[r_level];

      } else {

        if (start_run_rel_coords[r_level] ==
            domain.getGrid().getMinGridPoint(r_level)) {
          move_inverse.reset(r_level);
        } else if (end_run_rel_coords[r_level] ==
                   domain.getGrid().getMaxGridPoint(r_level)) {
          move_inverse.set(r_level);
        }

        if (move_inverse.test(r_level)) {

          // assert(runTypePos[r_level]>sl.getStartIndex(r_level,
          // startIndicesPos[s_level]));

          --runTypePos[r_level];

          end_run_rel_coords[r_level] = start_run_rel_coords[r_level] - 1;
          start_run_rel_coords[r_level] = sl.getRunStartCoord(
              r_level, startIndicesPos[s_level], runTypePos[r_level]);

          ++endRunAbsCoords[r_level];
          startRunAbsCoords[r_level] = endRunAbsCoords[r_level];

          if (start_run_rel_coords[r_level] ==
              domain.getGrid().getMinGridPoint(r_level)) {
            endRunAbsCoords[r_level] +=
                (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
          }

        } else {

          // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
          // startIndicesPos[s_level]+1)-1);

          ++runTypePos[r_level];

          start_run_rel_coords[r_level] = end_run_rel_coords[r_level] + 1;
          end_run_rel_coords[r_level] = sl.getRunEndCoord(
              r_level, startIndicesPos[s_level], runTypePos[r_level]);

          ++endRunAbsCoords[r_level];
          startRunAbsCoords[r_level] = endRunAbsCoords[r_level];

          if (end_run_rel_coords[r_level] ==
              domain.getGrid().getMaxGridPoint(r_level)) {
            endRunAbsCoords[r_level] +=
                (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
          }
        }
      }

      endRunAbsCoords[r_level] +=
          (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
      if (endRunAbsCoords[r_level] > domain.getGrid().getMaxGridPoint(r_level))
        endRunAbsCoords[r_level] = domain.getGrid().getMaxGridPoint(r_level);

      // assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
      // startIndicesPos[s_level]));       //
      // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
      // startIndicesPos[s_level]+1));

      // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

      // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // assert(startRunAbsCoords[r_level]<=endRunAbsCoords[r_level]);
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);

      if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
          end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
        // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
        // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
        // ||
        // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
      }

      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]);
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]);
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      return true;
    }

    return false;
  }

  bool go_previous_B() {

    // assert(s_level==r_level+1);
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxGridPoint(r_level));
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));

    const DomainSegmentType &sl = domain.domainSegments[sub];

    // assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
    // startIndicesPos[s_level]));     //
    // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
    // startIndicesPos[s_level]+1));

    if (startRunAbsCoords[r_level] !=
        domain.getGrid().getMinGridPoint(r_level)) {

      if (domain.getGrid().isBoundaryPeriodic(
              r_level)) { // periodic boundary conditions

        if (start_run_rel_coords[r_level] ==
            domain.getGrid().getMinGridPoint(r_level)) {
          runTypePos[r_level] =
              sl.getStartIndex(r_level, startIndicesPos[s_level] + 1) - 1;
          end_run_rel_coords[r_level] =
              domain.getGrid().getMaxGridPoint(r_level);
        } else {
          --runTypePos[r_level];
          end_run_rel_coords[r_level] = start_run_rel_coords[r_level] - 1;
        }

        start_run_rel_coords[r_level] = sl.getRunStartCoord(
            r_level, startIndicesPos[s_level], runTypePos[r_level]);

        --startRunAbsCoords[r_level];
        endRunAbsCoords[r_level] = startRunAbsCoords[r_level];

      } else {

        if (start_run_rel_coords[r_level] ==
            domain.getGrid().getMinGridPoint(r_level)) {
          move_inverse.set(r_level);
        } else if (end_run_rel_coords[r_level] ==
                   domain.getGrid().getMaxGridPoint(r_level)) {
          move_inverse.reset(r_level);
        }

        if (move_inverse.test(r_level)) {

          // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
          // startIndicesPos[s_level]+1)-1);

          ++runTypePos[r_level];

          start_run_rel_coords[r_level] = end_run_rel_coords[r_level] + 1;
          end_run_rel_coords[r_level] = sl.getRunEndCoord(
              r_level, startIndicesPos[s_level], runTypePos[r_level]);

          --startRunAbsCoords[r_level];
          endRunAbsCoords[r_level] = startRunAbsCoords[r_level];

          if (end_run_rel_coords[r_level] ==
              domain.getGrid().getMaxGridPoint(r_level)) {
            startRunAbsCoords[r_level] -=
                (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
          }

        } else {

          // assert(runTypePos[r_level]>sl.getStartIndex(r_level,
          // startIndicesPos[s_level]));

          --runTypePos[r_level];

          end_run_rel_coords[r_level] = start_run_rel_coords[r_level] - 1;
          start_run_rel_coords[r_level] = sl.getRunStartCoord(
              r_level, startIndicesPos[s_level], runTypePos[r_level]);

          --startRunAbsCoords[r_level];
          endRunAbsCoords[r_level] = startRunAbsCoords[r_level];

          if (start_run_rel_coords[r_level] ==
              domain.getGrid().getMinGridPoint(r_level)) {
            startRunAbsCoords[r_level] -=
                (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
          }
        }
      }

      startRunAbsCoords[r_level] -=
          (end_run_rel_coords[r_level] - start_run_rel_coords[r_level]);
      if (startRunAbsCoords[r_level] <
          domain.getGrid().getMinGridPoint(r_level))
        startRunAbsCoords[r_level] = domain.getGrid().getMinGridPoint(r_level);

      // assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
      // startIndicesPos[s_level]));
      // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
      // startIndicesPos[s_level]+1));

      // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

      // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // assert(startRunAbsCoords[r_level]<=endRunAbsCoords[r_level]);
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);

      if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
          end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
        // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
        // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
        // ||
        // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
      }

      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]);
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]);
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      return true;
    }

    return false;
  }

public:
  using DomainType = hrleDomain;

#ifdef NDEBUG
  void print() const {}
#else
  void print() const {
    BaseIterator<hrleDomain>::print();
    std::cout << "start_run_rel_coords: " << start_run_rel_coords << std::endl;
    std::cout << "end_run_rel_coords: " << end_run_rel_coords << std::endl;
    std::cout << "rel_coords: " << rel_coords << std::endl;
    std::cout << "offset: " << offset << std::endl;
    std::cout << "move_inverse: " << move_inverse << std::endl;
  }
#endif

  template <class V>
  SparseOffsetIterator(hrleDomain &passedDomain, const V &o,
                       bool reverse = false)
      : BaseIterator<hrleDomain>(passedDomain), offset(o) {
    if (reverse) {
      goToIndices(domain.getGrid().getMaxGridPoint());
    } else {
      goToIndices(domain.getGrid().getMinGridPoint());
    }
  }

  template <class V1, class V2>
  SparseOffsetIterator(hrleDomain &passedDomain, const V1 &o, const V2 &v)
      : BaseIterator<hrleDomain>(passedDomain), offset(o) {
    goToIndices(v);
  }

  SparseOffsetIterator &operator++() {
    while (true) {
      if (go_next_A())
        break;
      go_up_AB();
      if (go_next_B()) {
        go_down_BA(startRunAbsCoords[r_level]);
        break;
      }
      go_up_BA();
      // assert(it.r_level==it.s_level);
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

    int s = BaseIterator<hrleDomain>::getSegmentRun();
    if (s != sub) {
      goToIndices(s, absCoords);
    }

    return *this;
  }

  SparseOffsetIterator operator++(int) {
    SparseOffsetIterator temp(*this);
    ++(*this);
    return temp;
  }

  SparseOffsetIterator &operator--() {
    while (true) {
      if (go_previous_A())
        break;
      go_up_AB();
      if (go_previous_B()) {
        go_down_BA(endRunAbsCoords[r_level]);
        break;
      }
      go_up_BA();
      // assert(it.r_level==it.s_level);
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

    int s = BaseIterator<hrleDomain>::getSegmentRun();
    if (s != sub) {
      goToIndices(s, BaseIterator<hrleDomain>::getEndIndices());
      /*if (s<it.l.segmentation.size()) {
      assert(it.l.grid().decrementIndices(it.l.segmentation[s])==it.getEndIndices());
      } else {
      assert(it.l.grid().max_point_index()==it.getEndIndices());
      }*/
    }
    return *this;
  }

  SparseOffsetIterator operator--(int) {
    SparseOffsetIterator temp(*this);
    --(*this);
    return temp;
  }

  /// safe version of operator++(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool next() {
    if (BaseIterator<hrleDomain>::isFinished())
      return false;
    ++(*this);
    return true;
  }

  /// safe version of operator--(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool previous() {
    if (BaseIterator<hrleDomain>::isFinished())
      return false;
    --(*this);
    return true;
  }

  template <class V> void goToIndices(const V &v) {
    goToIndices(0, v); // TODO
    int s = BaseIterator<hrleDomain>::getSegmentRun();
    // std::cout << "got_to_indices, sub: " << s << std::endl;
    if (s != 0)
      goToIndices(s, v);
  }

  template <class V> void goToIndices(int subDomain, const V &v) {
    r_level = D;
    s_level = D;
    sub = subDomain;
    startIndicesPos[D] = 0;
    do {
      // assert(it.r_level==it.s_level);
      const IndexType c = v[r_level - 1];
      // assert(isDefined(it.startIndicesPos[it.s_level]));
      go_down_AB(c);
      // assert(c>=it.startRunAbsCoords[it.r_level]);
      // assert(it.endRunAbsCoords[it.r_level]>=c);
      go_down_BA(c);
    } while (r_level == s_level && s_level > 0);
    // assert(!it.isFinished());

    for (int h = 0; h < r_level; ++h) {
      absCoords[h] = domain.getGrid().getMinGridPoint(h);
      endAbsCoords[h] = domain.getGrid().getMaxGridPoint(h);
    }

    // TODO: initialize absCoords
  }

  Index<D> getOffset() const { return offset; }

  Index<D> getOffsetIndices() const {
    Index<D> indices = BaseIterator<hrleDomain>::getStartIndices();
    auto &grid = domain.getGrid();
    auto gridMin = grid.getMinGridPoint();
    auto gridMax = grid.getMaxGridPoint();
    for (unsigned i = 0; i < D; ++i) {
      indices[i] += offset[i];
      if (indices[i] > gridMax[i]) {
        const auto lengthDiff = indices[i] - gridMax[i];
        if (grid.isBoundaryReflective(i)) {
          indices[i] = gridMax[i] - lengthDiff;
        } else {
          indices[i] = gridMin[i] + (lengthDiff - 1);
        }
      } else if (indices[i] < gridMin[i]) {
        const auto lengthDiff = indices[i] - gridMin[i];
        if (grid.isBoundaryReflective(i)) {
          indices[i] = gridMin[i] - lengthDiff;
        } else {
          indices[i] = gridMax[i] + (lengthDiff + 1);
        }
      }
    }
    return indices;
  }
};

// typedef for const iterator
template <class hrleDomain>
using ConstSparseOffsetIterator = SparseOffsetIterator<const hrleDomain>;

} // namespace viennahrle

#endif // HRLE_OFFSET_RUNS_ITERATOR_HPP
