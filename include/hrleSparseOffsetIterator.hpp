#ifndef HRLE_OFFSET_RUNS_ITERATOR_HPP
#define HRLE_OFFSET_RUNS_ITERATOR_HPP

#include <bitset>
#include <vector>

#include "hrleBaseIterator.hpp"

/// The hrleSparseOffsetIterator iterates over the runs exactly like
/// hrleSparseIterator, but is set at indices given by the offset it was
/// initialised with. So if the offset was (0,0,1) and this iterator is
/// currently at the point (2,2,2), calling getValue() will return the value at
/// (2,2,2), but calling indices() will return (2,2,3). By comparing its indices
/// with those of an hrleSparseIterator and only incrementing when necessary,
/// this iterator can be used to keep track of the cartesian neighbour of an
/// hrleSparseIterator, taking into account boundary conditions.
template <class hrleDomain>
class hrleSparseOffsetIterator : public hrleBaseIterator<hrleDomain> {

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

  hrleVectorType<hrleIndexType, D> start_run_rel_coords;
  hrleVectorType<hrleIndexType, D> end_run_rel_coords;
  hrleVectorType<hrleIndexType, D> rel_coords;
  hrleVectorType<hrleIndexType, D> offset;
  std::bitset<D> move_inverse;

protected:
  void go_down_AB(hrleIndexType abs_c) { // X         //find right runTypePos

    // shfdhsfhdskjhgf assert(s_level==r_level);

    // shfdhsfhdskjhgf assert(domain.isDefined(startIndicesPos[s_level]));

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    const hrleSizeType &s = startIndicesPos[s_level];

    --r_level;

    int cycles = 0;
    hrleIndexType rel_c = domain.getGrid().globalIndex2LocalIndex(
        r_level, abs_c, offset[r_level], cycles);

    // shfdhsfhdskjhgf assert(s<sl.startIndices[r_level].size());
    // shfdhsfhdskjhgf assert(r_level>=0);

    hrleSizeType r = sl.startIndices[r_level][s];

    typename std::vector<hrleIndexType>::const_iterator start_breaks =
        sl.runBreaks[r_level].begin() + (r - s);
    typename std::vector<hrleIndexType>::const_iterator end_breaks =
        sl.runBreaks[r_level].begin() +
        (sl.getStartIndex(r_level, s + 1) - (s + 1));

    // shfdhsfhdskjhgf assert(start_breaks<=end_breaks);

    typename std::vector<hrleIndexType>::const_iterator pos_breaks =
        std::upper_bound(start_breaks, end_breaks, rel_c);

    r += hrleSizeType(pos_breaks - start_breaks);

    runTypePos[r_level] = r;

    if (pos_breaks == start_breaks) {
      start_run_rel_coords[r_level] = domain.getGrid().getMinGridPoint(r_level);
    } else {
      // shfdhsfhdskjhgf assert(pos_breaks>start_breaks);
      start_run_rel_coords[r_level] = *(pos_breaks - 1);
    }

    if (pos_breaks == end_breaks) {
      end_run_rel_coords[r_level] = domain.getGrid().getMaxGridPoint(r_level);
    } else {
      // shfdhsfhdskjhgf assert(pos_breaks<end_breaks);
      end_run_rel_coords[r_level] = (*pos_breaks) - 1;
    }

    // shfdhsfhdskjhgf
    // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
    // shfdhsfhdskjhgf assert(start_run_rel_coords[r_level]<=rel_c);
    // shfdhsfhdskjhgf assert(rel_c<=end_run_rel_coords[r_level]);
    // shfdhsfhdskjhgf
    // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

    // calculate run_absCoords

    if (domain.getGrid().isBoundaryPeriodic(r_level)) {

      const hrleIndexType &rel_s = start_run_rel_coords[r_level];
      const hrleIndexType &rel_e = end_run_rel_coords[r_level];

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

        hrleIndexType rel_s = start_run_rel_coords[r_level];
        hrleIndexType rel_e = end_run_rel_coords[r_level];

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

    // shfdhsfhdskjhgf
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
    // shfdhsfhdskjhgf assert(startRunAbsCoords[r_level]<=abs_c);
    // shfdhsfhdskjhgf assert(abs_c<=endRunAbsCoords[r_level]);
    // shfdhsfhdskjhgf
    // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);

    if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
        end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
      // shfdhsfhdskjhgf
      // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
      // shfdhsfhdskjhgf
      // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
      // ||
      // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
    }

    // shfdhsfhdskjhgf
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
    // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
    // offset[r_level])<=end_run_rel_coords[r_level]);

    // shfdhsfhdskjhgf
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
    // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
    // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
    // offset[r_level])<=end_run_rel_coords[r_level]);
  }

  void go_down_AB_first() { // find right runTypePos

    // shfdhsfhdskjhgf assert(s_level==r_level);

    const hrleIndexType c = domain.getGrid().getMinGridPoint(r_level - 1);
    // shfdhsfhdskjhgf assert(l.isDefined(startIndicesPos[s_level]));
    go_down_AB(c);

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

  void go_down_AB_last() { // find right runTypePos

    // shfdhsfhdskjhgf assert(s_level==r_level);

    const hrleIndexType c = domain.getGrid().getMaxGridPoint(r_level - 1);
    // shfdhsfhdskjhgf assert(l.isDefined(startIndicesPos[s_level]));
    go_down_AB(c);

    // shfdhsfhdskjhgf assert(s_level>=1);
    // shfdhsfhdskjhgf assert(s_level==r_level+1);
  }

  bool go_down_BA(hrleIndexType abs_c) { // X

    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf assert(s_level>=1);

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    startIndicesPos[s_level - 1] = sl.runTypes[r_level][runTypePos[r_level]];
    if (sl.isPtIdDefined(startIndicesPos[s_level - 1])) {
      --s_level;
      // shfdhsfhdskjhgf assert(domain.isDefined(startIndicesPos[s_level]));

      hrleIndexType rel_c = domain.getGrid().globalIndex2LocalIndex(
          r_level, abs_c, offset[r_level]);

      // shfdhsfhdskjhgf assert(rel_c>=start_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf assert(rel_c<=end_run_rel_coords[r_level]);

      startIndicesPos[s_level] += (rel_c - start_run_rel_coords[r_level]);
      rel_coords[s_level] = rel_c;
      absCoords[s_level] = abs_c;
      endAbsCoords[s_level] = abs_c;

      // shfdhsfhdskjhgf assert(domain.isDefined(startIndicesPos[s_level]));

      return true;
    } else {
      absCoords[r_level] = startRunAbsCoords[r_level];
      endAbsCoords[r_level] = endRunAbsCoords[r_level];

      // shfdhsfhdskjhgf assert(s_level==r_level+1);
      return false;
    }
  }

  bool go_next_A() { // X

    if (s_level == r_level) {

      // shfdhsfhdskjhgf
      // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
      // shfdhsfhdskjhgf
      // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

      // if
      // (startIndicesPos[s_level]<domain.runtypes[r_level][runTypePos[r_level]]+(end_run_coords[r_level]-start_run_coords[r_level]))
      // {

      if (absCoords[s_level] < endRunAbsCoords[r_level]) {
        ++absCoords[s_level];
        ++endAbsCoords[s_level];

        if (domain.getGrid().isBoundaryPeriodic(r_level)) {
          // shfdhsfhdskjhgf
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
            // shfdhsfhdskjhgf assert(startIndicesPos[s_level]>0);
            --startIndicesPos[s_level];
            --rel_coords[s_level];
          } else {
            ++startIndicesPos[s_level];
            ++rel_coords[s_level];
          }
        }

        // shfdhsfhdskjhgf
        // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
        // shfdhsfhdskjhgf
        // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
        // shfdhsfhdskjhgf
        // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
        // shfdhsfhdskjhgf
        // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

        return true;
      }
    }
    return false;
  }

  bool go_previous_A() {

    if (s_level == r_level) {

      // shfdhsfhdskjhgf
      // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
      // shfdhsfhdskjhgf
      // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

      if (absCoords[s_level] > startRunAbsCoords[r_level]) {
        --absCoords[s_level];
        --endAbsCoords[s_level];

        if (domain.getGrid().isBoundaryPeriodic(r_level)) {
          // shfdhsfhdskjhgf
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

        // shfdhsfhdskjhgf
        // assert(rel_coords[s_level]<=end_run_rel_coords[r_level]);
        // shfdhsfhdskjhgf
        // assert(rel_coords[s_level]>=start_run_rel_coords[r_level]);
        // shfdhsfhdskjhgf
        // assert(absCoords[s_level]<=endRunAbsCoords[r_level]);
        // shfdhsfhdskjhgf
        // assert(absCoords[s_level]>=startRunAbsCoords[r_level]);

        return true;
      }
    }
    return false;
  }

  bool go_next_B() {

    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxGridPoint(r_level));
    // shfdhsfhdskjhgf
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    // shfdhsfhdskjhgf assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
    // startIndicesPos[s_level])); shfdhsfhdskjhgf
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

          // shfdhsfhdskjhgf
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

          // shfdhsfhdskjhgf
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

      // shfdhsfhdskjhgf assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
      // startIndicesPos[s_level])); shfdhsfhdskjhgf
      // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
      // startIndicesPos[s_level]+1));

      // shfdhsfhdskjhgf
      // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // shfdhsfhdskjhgf
      // assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

      // shfdhsfhdskjhgf
      // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // shfdhsfhdskjhgf
      // assert(startRunAbsCoords[r_level]<=endRunAbsCoords[r_level]);
      // shfdhsfhdskjhgf
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);

      if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
          end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
        // shfdhsfhdskjhgf
        // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
        // shfdhsfhdskjhgf
        // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
        // ||
        // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
      }

      // shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      // shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      return true;
    } else {
      return false;
    }
  }

  bool go_previous_B() {

    // shfdhsfhdskjhgf assert(s_level==r_level+1);
    // shfdhsfhdskjhgf
    // assert(endRunAbsCoords[r_level]<=domain.getGrid().getMaxGridPoint(r_level));
    // shfdhsfhdskjhgf
    // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));

    const hrleDomainSegmentType &sl = domain.domainSegments[sub];

    // shfdhsfhdskjhgf assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
    // startIndicesPos[s_level])); shfdhsfhdskjhgf
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

          // shfdhsfhdskjhgf
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

          // shfdhsfhdskjhgf
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

      // shfdhsfhdskjhgf assert(runTypePos[r_level]>=sl.getStartIndex(r_level,
      // startIndicesPos[s_level])); shfdhsfhdskjhgf
      // assert(runTypePos[r_level]<sl.getStartIndex(r_level,
      // startIndicesPos[s_level]+1));

      // shfdhsfhdskjhgf
      // assert(start_run_rel_coords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // shfdhsfhdskjhgf
      // assert(start_run_rel_coords[r_level]<=end_run_rel_coords[r_level]);
      // shfdhsfhdskjhgf
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=end_run_rel_coords[r_level]);

      // shfdhsfhdskjhgf
      // assert(startRunAbsCoords[r_level]>=domain.getGrid().getMinGridPoint(r_level));
      // shfdhsfhdskjhgf
      // assert(startRunAbsCoords[r_level]<=endRunAbsCoords[r_level]);
      // shfdhsfhdskjhgf
      // assert(domain.getGrid().getMaxGridPoint(r_level)>=endRunAbsCoords[r_level]);

      if (endRunAbsCoords[r_level] - startRunAbsCoords[r_level] >
          end_run_rel_coords[r_level] - start_run_rel_coords[r_level]) {
        // shfdhsfhdskjhgf
        // assert(!(domain.getGrid().isBoundaryPeriodic(r_level)));
        // shfdhsfhdskjhgf
        // assert(end_run_rel_coords[r_level]==domain.getGrid().getMaxGridPoint(r_level)
        // ||
        // start_run_rel_coords[r_level]==domain.getGrid().getMinGridPoint(r_level));
      }

      // shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,startRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      // shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])>=start_run_rel_coords[r_level]); shfdhsfhdskjhgf
      // assert(domain.getGrid().globalIndex2LocalIndex(r_level,endRunAbsCoords[r_level],
      // offset[r_level])<=end_run_rel_coords[r_level]);

      return true;
    } else {
      return false;
    }
  }

public:
  using DomainType = hrleDomain;

  void print() const {
    hrleBaseIterator<hrleDomain>::print();
    std::cout << "start_run_rel_coords: " << start_run_rel_coords << std::endl;
    std::cout << "end_run_rel_coords: " << end_run_rel_coords << std::endl;
    std::cout << "rel_coords: " << rel_coords << std::endl;
    std::cout << "offset: " << offset << std::endl;
    std::cout << "move_inverse: " << move_inverse << std::endl;
  }

  template <class V>
  hrleSparseOffsetIterator(hrleDomain &passedDomain, const V &o,
                           bool reverse = false)
      : hrleBaseIterator<hrleDomain>(passedDomain), offset(o) {
    if (reverse) {
      goToIndices(domain.getGrid().getMaxGridPoint());
    } else {
      goToIndices(domain.getGrid().getMinGridPoint());
    }
  }

  template <class V1, class V2>
  hrleSparseOffsetIterator(hrleDomain &passedDomain, const V1 &o, const V2 &v)
      : hrleBaseIterator<hrleDomain>(passedDomain), offset(o) {
    goToIndices(v);
  }

  hrleSparseOffsetIterator<hrleDomain> &operator++() {
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
      const hrleIndexType c = domain.getGrid().getMinGridPoint(r_level - 1);
      // go_down_AB(c);
      go_down_AB_first();
      go_down_BA(c);
    }

    int s = hrleBaseIterator<hrleDomain>::getSegmentRun();
    if (s != sub) {
      goToIndices(s, absCoords);
    }

    return *this;
  }

  hrleSparseOffsetIterator<hrleDomain> operator++(int) {
    hrleSparseOffsetIterator<hrleDomain> temp(*this);
    ++(*this);
    return temp;
  }

  hrleSparseOffsetIterator<hrleDomain> &operator--() {
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
      const hrleIndexType c = domain.getGrid().getMaxGridPoint(r_level - 1);
      // go_down_AB(c);
      go_down_AB_last();
      go_down_BA(c);
    }

    int s = hrleBaseIterator<hrleDomain>::getSegmentRun();
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

  hrleSparseOffsetIterator<hrleDomain> operator--(int) {
    hrleSparseOffsetIterator<hrleDomain> temp(*this);
    --(*this);
    return temp;
  }

  /// safe version of operator++(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool next() {
    if (hrleBaseIterator<hrleDomain>::isFinished())
      return false;
    ++(*this);
    return true;
  }

  /// safe version of operator--(), checks first, whether iterator is finished
  /// or not before incrementing. Returns false if iterator was already
  /// finished
  bool previous() {
    if (hrleBaseIterator<hrleDomain>::isFinished())
      return false;
    --(*this);
    return true;
  }

  template <class V> void goToIndices(const V &v) {
    goToIndices(0, v); // TODO
    int s = hrleBaseIterator<hrleDomain>::getSegmentRun();
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
      // shfdhsfhdskjhgf assert(it.r_level==it.s_level);
      const hrleIndexType c = v[r_level - 1];

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

    // TODO initialize absCoords
  }

  const hrleVectorType<hrleIndexType, D> getOffset() const { return offset; }

  const hrleVectorType<hrleIndexType, D> getOffsetIndices() const {
    hrleVectorType<hrleIndexType, D> indices =
        hrleBaseIterator<hrleDomain>::getStartIndices();
    auto &grid = domain.getGrid();
    auto gridMin = grid.getMinGridPoint();
    auto gridMax = grid.getMaxGridPoint();
    for (unsigned i = 0; i < D; ++i) {
      indices[i] += offset[i];
      if (indices[i] > gridMax[i]) {
        if(grid.isBoundaryReflective(i)){
          indices[i] = gridMax[i] - 1;
        }else{
          indices[i] = gridMin[i] + (indices[i] - gridMax[i] - 1);
        }
      } else if (indices[i] < gridMin[i]) {
        if(grid.isBoundaryReflective(i)){
          indices[i] = gridMin[i] + 1;
        }else{
          indices[i] = gridMax[i] + (indices[i] - gridMin[i] + 1);
        }
      }
    }
    return indices;
  }
};

// typedef for const iterator
template <class hrleDomain>
using hrleConstSparseOffsetIterator =
    hrleSparseOffsetIterator<const hrleDomain>;

#endif // HRLE_OFFSET_RUNS_ITERATOR_HPP
