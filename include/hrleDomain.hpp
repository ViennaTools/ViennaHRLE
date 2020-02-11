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
#include <istream>
#include <ostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hrleAllocationType.hpp"
#include "hrleDomainSegmentArray.hpp"
#include "hrleSizeType.hpp"
#include "hrleSparseIterator.hpp"
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
  template <class> friend class hrleSparseIterator;
  template <class> friend class hrleSparseOffsetIterator;

  static char getBitSizeOfNumber(int number) {
    number = std::abs(number);
    int numBits = 0;
    while (number != 0) {
      number /= 2;
      ++numBits;
    }
    return numBits;
  }

  static char getByteSizeOfNumber(int number) {
    int numBits = getBitSizeOfNumber(number);
    return numBits / 8 + (numBits % 8) ? 1 : 0;
  }

  template <class S, class Type,
            typename std::enable_if<std::is_arithmetic<Type>::value,
                                    std::nullptr_t>::type = nullptr>
  static S &writeValue(S &stream, Type &value) {
    stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
    return stream;
  }

  template <class S, class Type,
            typename std::enable_if<!std::is_arithmetic<Type>::value,
                                    std::nullptr_t>::type = nullptr>
  static S &writeValue(S &stream, Type &value) {
    value->serialize(stream);
    return stream;
  }

  template <class S, class Type,
            typename std::enable_if<std::is_arithmetic<Type>::value,
                                    std::nullptr_t>::type = nullptr>
  static S &readValue(S &stream, Type &value) {
    stream.read(reinterpret_cast<char *>(&value), sizeof(value));
    return stream;
  }

  template <class S, class Type,
            typename std::enable_if<!std::is_arithmetic<Type>::value,
                                    std::nullptr_t>::type = nullptr>
  static S &readValue(S &stream, Type &value) {
    value->deserialize(stream);
    return stream;
  }

public:
  // CONSTRUCTORS
  hrleDomain(){};

  // create empty level set with one undefined run
  hrleDomain(hrleGrid<D> &g, hrleSizeType runType = hrleRunTypeValues::UNDEF_PT)
      : grid(&g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
  };

  hrleDomain(hrleGrid<D> *g, hrleSizeType runType = hrleRunTypeValues::UNDEF_PT)
      : grid(g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
  };

  // create empty level set with one undefined value
  hrleDomain(hrleGrid<D> &g, T value) : grid(&g) {
    initialize();
    domainSegments[0].insertNextUndefinedPoint(grid->getMinIndex(), value);
  };

  hrleDomain(hrleGrid<D> *g, T value) : grid(g) {
    initialize();
    domainSegments[0].insertNextUndefinedPoint(grid->getMinIndex(), value);
  };

  void deepCopy(const hrleDomain<T, D> &passedDomain) {
    assert(grid == &(passedDomain.getGrid()));
    pointIdOffsets = passedDomain.pointIdOffsets;
    segmentation = passedDomain.segmentation;
    domainSegments.deepCopy(grid, passedDomain.domainSegments);
  }

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

  hrleGrid<D> &getGrid() { return *grid; }

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

  unsigned getNumberOfRuns(int segmentId, int dimension) const {
    return domainSegments[segmentId].getNumberOfRuns(dimension);
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
        bounds[2 * i] = grid->getMinGridPoint(i);
      }

      if (grid->getBoundaryConditions(i) == hrleGrid<D>::INFINITE_BOUNDARY ||
          grid->getBoundaryConditions(i) ==
              hrleGrid<D>::POS_INFINITE_BOUNDARY) {
        bounds[2 * i + 1] = getMaxRunBreak(i);
      } else {
        bounds[2 * i + 1] = grid->getMaxGridPoint(i);
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
      s.insertNextUndefinedRunType(segmentation.back(), grid->getMaxGridPoint(),
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

      hrleVectorType<hrleIndexType, D> startOfSegment =
          (p == 0) ? grid->getMinIndex() : newDomain.segmentation[p - 1];

      hrleVectorType<hrleIndexType, D> endOfSegment =
          (p != static_cast<int>(newDomain.segmentation.size()))
              ? newDomain.segmentation[p]
              : grid->getMaxGridPoint();

      for (hrleConstSparseIterator<hrleDomain> it(*this, startOfSegment);
           it.getStartIndices() < endOfSegment; it.next()) {
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
  void desegment() {
    if (domainSegments.getNumberOfSegments() > 1) {
      hrleDomain newDomain(grid);
      newDomain.initialize(); // initialize with only one segment

      hrleDomainSegmentType &s = newDomain.domainSegments[0];
      for (hrleConstSparseIterator<hrleDomain> it(*this); !it.isFinished();
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

  /*
  Serialized hrle structure
  *********************************************************************
  *    H-RLE Block Header: 15 Bytes    *
  **************************************
  * NOTE: The following is repeated for each dimension
  * 1 Byte    This byte contains the number of bytes
              used for each start index
  * 1 Byte    This byte contains the number of bits(!!!!!!)
              used for each run type
  * 1 Byte    This byte contains the number of bytes
              used for each runbreak
  * 4 Bytes   Number of saved Start Indices
  * 4 Bytes   Number of saved Runtypes
  * 4 Bytes   Number of saved Runbreaks
  *********************************************************************
  *    H-RLE Block Data    *
  **************************
  * Start Indices   using adaptive number of bytes(delta encoded)
  * Runtypes        using n bits per runtype
  * Indices of runtypes  using adaptive number of bytes(delta encoded)
  * Runbreaks       using adaptive number of bytes
  * NOTE: END
  *********************************************************************
  *    Values Header: 8 Bytes    *
  ***********************************
  * 4 Bytes  Number of defined values
  * 4 Bytes  Number of undefined values
  *********************************************************************
  *    Values Data    *
  ************************
  * Defined Values - using sizeOf Bytes per value or
                     serialize/deserialize of underlying data structure
  * Undefined Values - same as Defined Values
  *********************************************************************
  */
  /// Serialize the stored hrle structure.
  /// If a non-arithmetic type is used, it needs to implement
  /// the function "std::ostream &serialize(std::ostream)".
  /// If a parallelized hrle structure is used, it is segmented
  /// anew since only a serial (non-parallelized) hrle structure
  /// can be serialized.
  std::ostream &serialize(std::ostream &stream) {
    // write identifier
    stream << "hrleDomain";

    bool structureIsSerial = true;
    if (getNumberOfSegments() > 1) {
      desegment();
      structureIsSerial = false;
    }

    // get domain bounds
    hrleIndexType bounds[2 * D];
    getDomainBounds(bounds);

    // START INDICES, RUN TYPES, RUN BREAKS FOR EACH DIMENSION
    for (int dim = D - 1; dim >= 0; --dim) {
      // get start indices, runbreaks and runtypes
      const std::vector<hrleSizeType> &startIndices =
          domainSegments[0].startIndices[dim];
      const std::vector<hrleSizeType> &runTypes =
          domainSegments[0].runTypes[dim];
      const std::vector<hrleIndexType> &runBreaks =
          domainSegments[0].runBreaks[dim];

      const char bytesPerIndex =
          getByteSizeOfNumber(bounds[2 * dim + 1] - bounds[2 * dim]);

      char bitsPerRunType = getBitSizeOfNumber(getNumberOfUndefinedValues());
      if (bitsPerRunType == 3) {
        ++bitsPerRunType;
      } else if (bitsPerRunType > 4 && bitsPerRunType < 8) {
        bitsPerRunType = 8;
      } else if (bitsPerRunType > 8 && (bitsPerRunType % 8)) {
        bitsPerRunType += 8 - bitsPerRunType % 8;
      }

      const char bytesPerRunBreak =
          std::max(getBitSizeOfNumber(bounds[2 * dim]),
                   getBitSizeOfNumber(bounds[2 * dim + 1])) /
              8 +
          1;

      // HRLE BLOCK HEADER
      stream.write((char *)&bytesPerIndex, 1);
      stream.write((char *)&bitsPerRunType, 1);
      stream.write((char *)&bytesPerRunBreak, 1);
      {
        uint32_t numberOfValues = uint32_t(startIndices.size());
        stream.write((char *)&numberOfValues, 4);
        numberOfValues = uint32_t(runTypes.size());
        stream.write((char *)&numberOfValues, 4);
        numberOfValues = uint32_t(runBreaks.size());
        stream.write((char *)&numberOfValues, 4);
      }

      // uint32_t values_written = 0;
      // Write start indices; only save the difference to the next start index
      // (delta encoding)

      // First index is always 0, no need to write explicitly
      for (unsigned int i = 1; i < startIndices.size(); i++) {
        unsigned long diff = startIndices[i] - startIndices[i - 1];
        stream.write((char *)&diff, bytesPerIndex);
        // values_written++;
      }

      // write all runtypes to the file, skipping all segments and indices
      int count = 8 / bitsPerRunType - 1;
      unsigned char byte = 0;
      std::vector<hrleSizeType>
          definedRunIndices; // store all indices of defined runtypes

      // each runType needs at least one byte
      if (bitsPerRunType > 4) {
        for (typename std::vector<hrleSizeType>::const_iterator it =
                 runTypes.begin();
             it != runTypes.end(); ++it) {
          hrleSizeType PtId = 0;
          // if undefined point, need to shift id
          if (*it >= hrleRunTypeValues::UNDEF_PT)
            PtId = (*it) - hrleRunTypeValues::UNDEF_PT + 1;
          else
            definedRunIndices.push_back(*it);

          stream.write((char *)&PtId, (bitsPerRunType - 1) / 8 + 1);
        }
      } else { // can fit more than one value in a byte
        for (typename std::vector<hrleSizeType>::const_iterator it =
                 runTypes.begin();
             it != runTypes.end(); ++it) {
          hrleSizeType PtId = 0;
          if (*it >= hrleRunTypeValues::UNDEF_PT)
            PtId = (*it) - hrleRunTypeValues::UNDEF_PT + 1;
          else
            definedRunIndices.push_back(*it);

          byte |= (PtId << (count * bitsPerRunType));
          --count;

          if (count < 0) { // push byte to stream and start again
            stream << byte;
            count = 8 / bitsPerRunType - 1;
            byte = 0;
          }
        }
        // if last byte is not completely filled, just push it
        if (count >= 0 && count != 8 / bitsPerRunType - 1) {
          stream << byte;
        }
      }

      // Write indices of defined runtypes; only save the difference to the next
      // defined runtype write the first runtype(always 0) explicitly; makes
      // reading easier
      stream.write((char *)&definedRunIndices[0], bytesPerIndex);
      for (unsigned int i = 0; i < definedRunIndices.size() - 1; i++) {
        unsigned long diff = definedRunIndices[i + 1] - definedRunIndices[i];
        stream.write((char *)&diff, bytesPerIndex);
      }

      // Write runbreaks
      for (typename std::vector<hrleIndexType>::const_iterator it =
               runBreaks.begin();
           it != runBreaks.end(); ++it) {
        stream.write((char *)&(*it), bytesPerRunBreak);
      }
    }

    // DEFINED VALUES, UNDEFINED VALUES
    // IMPORTANT NOTE:
    // If a non-arithmetic type is used, it needs to implement
    // the function "std::ostream &serialize(std::ostream)"
    // HEADER
    {
      uint32_t numberOfDefinedValues = getNumberOfPoints();
      uint32_t numberOfUndefinedValues = getNumberOfUndefinedValues();
      stream.write((char *)&numberOfDefinedValues, 4);
      stream.write((char *)&numberOfUndefinedValues, 4);
    }
    // DATA
    {
      std::vector<hrleValueType> &definedValues =
          domainSegments[0].definedValues;
      for (typename std::vector<hrleValueType>::const_iterator it =
               definedValues.begin();
           it != definedValues.end(); ++it) {
        writeValue(stream, *it);
        // stream.write(reinterpret_cast<char*>*it) //<< *it;
      }

      std::vector<hrleValueType> &undefinedValues =
          domainSegments[0].undefinedValues;
      for (typename std::vector<hrleValueType>::const_iterator it =
               undefinedValues.begin();
           it != undefinedValues.end(); ++it) {
        writeValue(stream, *it);
        // stream << *it;
      }
    }

    // if the hrleDomain was segmented before, segment it again
    if (!structureIsSerial) {
      segment();
    }

    return stream;
  }

  std::istream &deserialize(std::istream &stream) {
    // check identifier
    char identifier[10];
    stream.read(identifier, 10);
    if (std::string(identifier).compare(0, 10, "hrleDomain")) {
      std::cout
          << "Reading hrleDomain from stream failed. Header could not be found."
          << std::endl;
      return stream;
    }

    // READ HRLE PROPERTIES
    for (int dim = D - 1; dim >= 0; --dim) {
      // get the start indices, runtypes and runbreaks vectors
      std::vector<hrleSizeType> &startIndices =
          domainSegments[0].startIndices[dim];
      std::vector<hrleSizeType> &runTypes = domainSegments[0].runTypes[dim];
      std::vector<hrleIndexType> &runBreaks = domainSegments[0].runBreaks[dim];

      uint32_t numberOfStartIndices, numberOfRunTypes, numberOfRunBreaks;
      char bytesPerIndex, bitsPerRunType, bytesPerRunBreak;
      // reading in the 15 byte H-RLE header
      stream.read(&bytesPerIndex, 1);
      stream.read(&bitsPerRunType, 1);
      stream.read(&bytesPerRunBreak, 1);
      stream.read((char *)&numberOfStartIndices, 4);
      stream.read((char *)&numberOfRunTypes, 4);
      stream.read((char *)&numberOfRunBreaks, 4);

      // READ START INDICES
      {
        startIndices.clear();
        unsigned long long sum = 0;
        // push the 0, it was not written to the file
        startIndices.push_back(0);
        for (unsigned i = 0; i < numberOfStartIndices - 1; ++i) {
          unsigned long long current = 0;
          stream.read((char *)&current, bytesPerIndex);
          sum += current;
          startIndices.push_back(hrleSizeType(sum));
        }
      }

      // READ RUN TYPES
      {
        runTypes.clear();
        unsigned numberOfValuesPerByte = 8 / bitsPerRunType;
        unsigned numberOfBytes =
            (numberOfRunTypes - 1) / numberOfValuesPerByte + 1;
        // Read defined run IDs with second file stream
        // std::ifstream runIdFin(filePath);
        // runIdFin.seekg((long)stream.tellg() + (long)numberOfBytes);

        // unsigned long long definedId = 0;
        runTypes.resize(numberOfRunTypes);
        if (bitsPerRunType > 4) {
          for (unsigned i = 0; i < numberOfRunTypes; ++i) {
            hrleSizeType current = 0;
            stream.read((char *)&current, (bitsPerRunType - 1) / 8 + 1);
            if (current == 0) { // defined run
              // for a defined run, store 0, so it can be replaced later
              runTypes[i] = 0;
            } else { // undefined run
              runTypes[i] = current - 1 + hrleRunTypeValues::UNDEF_PT;
            }
          }
        } else { // if there are several values in each byte
          unsigned readValues = 0;
          unsigned long long bitMask = ((1 << bitsPerRunType) - 1)
                                       << (8 - bitsPerRunType);
          for (unsigned i = 0; i < numberOfBytes; ++i) {
            char byte;
            stream.read(&byte, 1);
            for (unsigned j = 0; j < numberOfValuesPerByte; ++j) {
              if (readValues == numberOfRunTypes)
                break;
              hrleSizeType current = (byte & bitMask) >> (8 - bitsPerRunType);
              byte <<= bitsPerRunType;
              if (current == 0) { // defined run
                runTypes[i * numberOfValuesPerByte + j] = 0;
              } else { // undefined run
                runTypes[i * numberOfValuesPerByte + j] =
                    current - 1 + hrleRunTypeValues::UNDEF_PT;
              }
              ++readValues;
            }
          }
        }
        assert(runTypes.size() == numberOfRunTypes);

        // READ DEFINED RUNS
        // they are delta encoded so, add them onto definedId
        {
          unsigned long long definedId = 0;
          unsigned long long counter = 0;
          while (counter < runTypes.size()) {
            if (runTypes[counter] == 0) {
              unsigned long long relativeId = 0;
              stream.read(reinterpret_cast<char *>(&relativeId), bytesPerIndex);
              definedId += relativeId;
              runTypes[counter] = definedId;
            }
            ++counter;
          }
        }

        // stream.seekg(runIdFin.tellg());
        // runIdFin.close();
      }

      // READ RUNBREAKS
      {
        runBreaks.clear();
        // bitmask sets all higher bytes to FF for negative numbers
        long long bitMask = 0;
        --bitMask <<= bytesPerRunBreak * 8;
        for (unsigned i = 0; i < numberOfRunBreaks; ++i) {
          long long runBreak = 0;
          stream.read((char *)&runBreak, bytesPerRunBreak);
          if (runBreak >> (bytesPerRunBreak * 8 - 1))
            runBreak |= bitMask;
          runBreaks.push_back(hrleSizeType(runBreak));
        }
      }
    }

    // READ DEFINED VALUES
    {
      // HEADER
      unsigned numberOfDefinedValues = 0, numberOfUndefinedValues = 0;
      stream.read((char *)&numberOfDefinedValues, 4);
      stream.read((char *)&numberOfUndefinedValues, 4);

      domainSegments[0].definedValues.clear();
      // DEFINED VALUES
      for (unsigned i = 0; i < numberOfDefinedValues; ++i) {
        hrleValueType definedValue;
        readValue(stream, definedValue);
        domainSegments[0].definedValues.push_back(definedValue);
      }

      domainSegments[0].undefinedValues.clear();
      // UNDEFINED VALUES
      for (unsigned i = 0; i < numberOfUndefinedValues; ++i) {
        hrleValueType undefinedValue;
        readValue(stream, undefinedValue);
        domainSegments[0].undefinedValues.push_back(undefinedValue);
      }
    }

    finalize();
    segment();

    return stream;
  }
};

#endif // HRLE_DOMAIN_HPP_
