#ifndef HRLE_DOMAIN_HPP_
#define HRLE_DOMAIN_HPP_

/* =========================================================================
Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

-----------------
ViennaTS - The Vienna Topography Simulator
-----------------

Contact:         viennats@iue.tuwien.ac.at

License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#include <algorithm>
#include <bitset>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hrleAllocationType.hpp"
#include "hrleDomainSegmentArray.hpp"
#include "hrleRunsIterator.hpp"
#include "hrleSizeType.hpp"
#include "hrleVectorType.hpp"

static constexpr unsigned int lvlset_omp_max_num_threads = 1000;

template <class T = double, int D = 3> class hrleDomain {
public:
  // TYPEDEFS
  typedef hrleDomainSegment<T, D> hrleDomainSegmentType;
  typedef typename hrleDomainSegmentType::hrleValueType hrleValueType;

  typedef hrleVectorType<hrleIndexType, D>
      hrleIndexPoint; // 3d index vector type
  typedef std::vector<hrleIndexPoint> hrleIndexPoints;

  // CONSTANTS
  static constexpr int dimension = D;

private:
  hrleDomain(hrleDomain &){};

  // Private member variables
  hrleGrid<D> *grid; // Grid stores the information about the grid, on
  // which the HRLE structure is defined

  hrleIndexPoints segmentation;
  hrleDomainSegmentArray<T, D> domainSegments;

  // std::vector<hrleSizeType> activePointIdOffsets;
  std::vector<hrleSizeType> pointIdOffsets;

  // Friend classes
  // Iterators need private access to get the values
  template <class> friend class hrleBaseIterator;
  template <class> friend class hrleRunsIterator;
  template <class> friend class hrleOffsetRunsIterator;
  // reader and writer need private acces for runtypes/runbreaks/etc
  template <class> friend class hrleDomainWriter;
  template <class> friend class hrleDomainReader;

public:
  // CONSTRUCTORS
  hrleDomain(){};

  // create empty level set with one undefined run
  hrleDomain(hrleGrid<D> &g, T runType = hrleRunTypeValues::UNDEF_PT)
      : grid(&g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
    // finalize();
  };

  hrleDomain(hrleGrid<D> *g, T runType = hrleRunTypeValues::UNDEF_PT)
      : grid(g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
    // finalize();
  };

  void deepCopy(hrleGrid<D> *passedGrid, const hrleDomain<T, D> &passedDomain) {
    deepCopy(*passedGrid, passedDomain);
  }

  void deepCopy(hrleGrid<D> &passedGrid, const hrleDomain<T, D> &passedDomain) {
    grid = &passedGrid;
    pointIdOffsets = passedDomain.pointIdOffsets;
    segmentation = passedDomain.segmentation;
    domainSegments.deepCopy(grid, passedDomain.domainSegments);
  }

  void shallowCopy(const hrleDomain<T, D> &passedDomain) {
    grid = passedDomain.grid;
    pointIdOffsets = passedDomain.pointIdOffsets;
    domainSegments.shallowCopy(passedDomain.domainSegments);
    segmentation = passedDomain.segmentation;
  }

  // Inline member functions
  void print(std::ostream &out = std::cout) const {
    for (hrleSizeType i = 0; i != domainSegments.getNumberOfSegments(); ++i) {
      domainSegments[i].print(out);
      out << std::endl;
    }
  }

  const hrleGrid<D> &getGrid() const { return *grid; }

  hrleDomainSegmentType &getDomainSegment(unsigned i) {
    return domainSegments[i];
  }

  hrleSizeType getPointIdOffset(unsigned i) const { return pointIdOffsets[i]; }

  unsigned getNumberOfPoints() const {
    unsigned totalPoints = 0;
    for (unsigned i = 0; i < domainSegments.getNumberOfSegments(); ++i) {
      totalPoints += domainSegments[i].getNumberOfPoints();
    }
    return totalPoints;
  }

  unsigned getNumberOfUndefinedValues() const {
    unsigned total = 0;
    for (unsigned i = 0; i < domainSegments.getNumberOfSegments(); ++i) {
      total += domainSegments[i].getNumberOfUndefinedValues();
    }
    return total;
  }

  unsigned getNumberOfSegments() const {
    return unsigned(domainSegments.getNumberOfSegments());
  }

  const hrleIndexPoints &getSegmentation() const { return segmentation; }

  hrleIndexType getMinRunBreak(int dim) const {
    hrleIndexType minBreak = domainSegments[0].getMinRunBreak(dim);
    for (unsigned i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      minBreak = std::min(minBreak, domainSegments[i].getMinRunBreak(dim));
    }
    return minBreak;
  }

  hrleIndexType getMaxRunBreak(int dim) const {
    hrleIndexType maxBreak = domainSegments[0].getMaxRunBreak(dim);
    for (unsigned i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      maxBreak = std::max(maxBreak, domainSegments[i].getMaxRunBreak(dim));
    }
    return maxBreak;
  }

  hrleVectorType<hrleIndexType, D> getMinRunBreak() const {
    hrleVectorType<hrleIndexType, D> minBreak;
    for (unsigned i = 0; i < D; ++i) {
      minBreak[i] = getMinRunBreak(i);
    }
    return minBreak;
  }

  hrleVectorType<hrleIndexType, D> getMaxRunBreak() const {
    hrleVectorType<hrleIndexType, D> maxBreak;
    for (unsigned i = 0; i < D; ++i) {
      maxBreak[i] = getMaxRunBreak(i);
    }
    return maxBreak;
  }

  void getDomainBounds(hrleIndexType *bounds) {
    for (unsigned i = 0; i < D; ++i) {
      if (grid->getBoundaryConditions(i) == hrleGrid<D>::INFINITE_BOUNDARY ||
          grid->getBoundaryConditions(i) ==
              hrleGrid<D>::NEG_INFINITE_BOUNDARY) {
        bounds[2 * i] = getMinRunBreak(i);
      } else {
        bounds[2 * i] = grid->getMinIndex(i);
      }

      if (grid->getBoundaryConditions(i) == hrleGrid<D>::INFINITE_BOUNDARY ||
          grid->getBoundaryConditions(i) ==
              hrleGrid<D>::POS_INFINITE_BOUNDARY) {
        bounds[2 * i + 1] = getMaxRunBreak(i);
      } else {
        bounds[2 * i + 1] = grid->getMaxIndex(i);
      }
    }
  }

  template <class V>
  void insertNextDefinedPoint(int sub, const V &point, T value) {
    // this function adds a new defined grid point to the level set function
    // the indices of the new grid point "indices" must be greater
    // (lexiographically greater) than the indices of the current "last" grid
    // point
    domainSegments[sub].insertNextDefinedPoint(point, value);
  }

  template <class V>
  void insertNextUndefinedRunType(int sub, const V &point, hrleSizeType rt) {
    // this function sets the sign of an undefined run starting at indices
    // "point" the indices of the new grid point "indices" must be greater
    // (lexiographically greater) than the indices of the current "last" grid
    // point
    domainSegments[sub].insertNextUndefinedRunType(point, rt);
  }

  template <class V>
  void insertNextUndefinedPoint(int sub, const V &point, hrleValueType value) {
    // this function sets the sign of an undefined run starting at indices
    // "point" the indices of the new grid point "indices" must be greater
    // (lexiographically greater) than the indices of the current "last" grid
    // point
    domainSegments[sub].insertNextUndefinedPoint(point, value);
  }

  // create new segmentation and domainSegments
  void initialize(const hrleIndexPoints &p = hrleIndexPoints(),
                  const hrleAllocationType<hrleSizeType, D> &a =
                      hrleAllocationType<hrleSizeType, D>()) {

    segmentation = p;

    domainSegments.initialize(segmentation.size() + 1, *grid,
                              a); // clear old segments and create new ones

    for (hrleSizeType i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      hrleDomainSegment<T, D> &s = domainSegments[i];

      s.insertNextUndefinedRunType(grid->getMinIndex(),
                                   grid->decrementIndices(segmentation[0]),
                                   hrleRunTypeValues::SEGMENT_PT);

      for (hrleSizeType j = 1; j < i; ++j) {
        s.insertNextUndefinedRunType(segmentation[j - 1],
                                     grid->decrementIndices(segmentation[j]),
                                     hrleRunTypeValues::SEGMENT_PT + j);
      }
    }
  }

  // distribute data points evenly across hrleDomainSegments and add SEGMENT_PT
  // as boundary markers
  void finalize() {

    for (hrleSizeType i = 0; i + 1 < domainSegments.getNumberOfSegments();
         ++i) {

      hrleDomainSegment<T, D> &s = domainSegments[i];

      for (hrleSizeType j = i + 1; j < segmentation.size(); ++j) {
        s.insertNextUndefinedRunType(segmentation[j - 1],
                                     grid->decrementIndices(segmentation[j]),
                                     hrleRunTypeValues::SEGMENT_PT + j);
      }
      s.insertNextUndefinedRunType(segmentation.back(), grid->getMaxIndex(),
                                   hrleRunTypeValues::SEGMENT_PT +
                                       hrleSizeType(segmentation.size()));
    }

    // TODO: for now do not save pointId offsets as they can easily be found
    // from cumulating domainSegments[i].getNumberOfPoints()
    // calculate id-offsets
    pointIdOffsets.clear();
    // active_pointIdOffsets.clear();
    pointIdOffsets.push_back(0);
    // activePointId_offsets.push_back(0);
    for (hrleSizeType i = 0; i < domainSegments.getNumberOfSegments() - 1;
         ++i) {
      pointIdOffsets.push_back(pointIdOffsets.back() +
                               domainSegments[i].getNumberOfPoints());
      //     activePointId_offsets.push_back(activePointId_offsets.back()+domainSegments[i].num_active_pts());
    }

    // assert(activePointId_offsets.size()==hrleDomainSegment.size());
    assert(pointIdOffsets.size() == domainSegments.getNumberOfSegments());
  }

  /// Converts a pointId(given by lexicographical order of points) to a spatial
  /// coordinate
  hrleIndexPoint ptIdToCoordinate(hrleSizeType pt) const {
    hrleIndexPoint pointCoords;

    // find domainSegment of point
    const int segment = int(
        std::upper_bound(pointIdOffsets.begin() + 1, pointIdOffsets.end(), pt) -
        (pointIdOffsets.begin() + 1));

    pt -= pointIdOffsets[segment]; // local point id

    const hrleDomainSegmentType &s = domainSegments[segment];

    // find right PointID by bisection
    for (int g = 0; g < D; ++g) {
      hrleSizeType min = 0;
      hrleSizeType max = hrleSizeType(s.runTypes[g].size()) - 1;

      while (!s.isPtIdDefined(s.runTypes[g][min]))
        ++min;
      while (!s.isPtIdDefined(s.runTypes[g][max]))
        --max;

      while (min != max) {
        hrleSizeType mid = (max + min + 1) / 2;
        while (!s.isPtIdDefined(s.runTypes[g][mid]))
          ++mid;

        if (s.runTypes[g][mid] <= pt) {
          min = mid;
        } else {
          max = (max + min - 1) / 2;
          while (!s.isPtIdDefined(s.runTypes[g][max]))
            --max;
        }
      }

      pointCoords[g] = pt - s.runTypes[g][min];

      pt = hrleSizeType(std::upper_bound(s.startIndices[g].begin() + 1,
                                         s.startIndices[g].end(), min) -
                        (s.startIndices[g].begin() + 1));

      pointCoords[g] += s.getRunStartCoord(g, pt, min);
    }

    return pointCoords;
  }

  /// finds the ideal coordinates at which to break between domainSegments
  /// for balanced distribution of points
  hrleIndexPoints getNewSegmentation() const {
    hrleIndexPoints tempSegmentation;

    int n = 1;
#ifdef _OPENMP
    n = omp_get_max_threads();
#endif
    hrleSizeType n_pts = getNumberOfPoints(); // number of defined grid points
    hrleSizeType sum = 0;

    for (unsigned int g = 0; g < static_cast<unsigned int>(n - 1); ++g) {
      sum += n_pts / n + ((n_pts % n) > g);
      // TODO: is that if really necessary?
      if (sum != n_pts)
        tempSegmentation.push_back(ptIdToCoordinate(sum));
    }

    return tempSegmentation;
  }

  /// allocation_type allocates the requried sizes to num_values and num_runs
  /// for all the domainSegments members num_values[0] is to contain level set
  /// values, num_values[i] contains the start indices at the i-th dimension
  /// num_runs[i] is to contain the run types at the i-th dimension
  hrleAllocationType<hrleSizeType, D> getAllocation() const {
    hrleAllocationType<hrleSizeType, D> a;
    a.num_values = hrleVectorType<hrleIndexType, D>(0);
    a.num_runs = hrleVectorType<hrleIndexType, D>(0);
    for (unsigned i = 0; i < domainSegments.getNumberOfSegments(); ++i) {
      hrleAllocationType<hrleSizeType, D> b = domainSegments[i].getAllocation();
      a.num_values = Max(a.num_values, b.num_values);
      a.num_runs = Max(a.num_runs, b.num_runs);
    }

    return (a * domainSegments.getNumberOfSegments());
  }

  /// distribute points evenly across domainSegments, so that they can be
  /// iterated over by separate threads efficiently
  void segment() {
    hrleDomain newDomain(grid);
    newDomain.initialize(getNewSegmentation(), getAllocation());

#pragma omp parallel num_threads(                                              \
    int(newDomain.domainSegments.getNumberOfSegments()))
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      hrleDomainSegmentType &s = newDomain.domainSegments[p];
      hrleConstRunsIterator<hrleDomain> it(
          *this,
          (p == 0) ? grid->getMinIndex() : newDomain.segmentation[p - 1]);

      hrleVectorType<hrleIndexType, D> endOfSegment =
          (p != static_cast<int>(newDomain.segmentation.size()))
              ? newDomain.segmentation[p]
              : grid->getMaxIndex();

      for (; it.getStartIndices() < endOfSegment; it.next()) {
        if (it.isDefined()) {
          s.insertNextDefinedPoint(it.getStartIndices(), it.getValue());
        } else {
          s.insertNextUndefinedPoint(it.getStartIndices(), it.getValue());
        }
      }
    }

    newDomain.finalize();
    deepCopy(*grid, newDomain);
  }

  /// opposite of "segment()" puts all data back into a single
  /// HRLE structure
  void serialize() {
    if (domainSegments.getNumberOfSegments() > 1) {
      hrleDomain newDomain(grid);
      newDomain.initialize(); // initialize with only one segment

      hrleDomainSegmentType &s = newDomain.domainSegments[0];
      for (hrleConstRunsIterator<hrleDomain> it(*this); !it.isFinished();
           ++it) {
        if (it.isDefined()) {
          s.insertNextDefinedPoint(it.getStartIndices(), it.getDefinedValue());
        } else {
          s.insertNextUndefinedPoint(it.getStartIndices(), it.getValue());
        }
      }

      newDomain.finalize();
      deepCopy(*grid, newDomain);
    }
  }
};

#endif // HRLE_DOMAIN_HPP_
