#ifndef HRLE_DOMAIN_HPP_
#define HRLE_DOMAIN_HPP_

/* =========================================================================
Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

-----------------
ViennaHRLE
-----------------

Contact:         viennatools@iue.tuwien.ac.at

License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */

#include <algorithm>
#include <bitset>
#include <limits>
#include <string>
#include <type_traits>
#include <vector>

#include <iostream>
#include <istream>
#include <ostream>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "hrleAllocationType.hpp"
#include "hrleDomainSegmentArray.hpp"
#include "hrleSparseIterator.hpp"
#include "hrleTypes.hpp"
#include "hrleUtil.hpp"

namespace viennahrle {
using namespace viennacore;
template <class T = double, int D = 3> class Domain {
public:
  // TYPEDEFS
  typedef DomainSegment<T, D> DomainSegmentType;
  typedef typename DomainSegmentType::ValueType ValueType;
  typedef std::vector<Index<D>> IndexPoints;

  // CONSTANTS
  static constexpr int dimension = D;

private:
  // Private member variables
  Grid<D> *grid; // Grid stores the information about the grid, on
  // which the HRLE structure is defined

  IndexPoints segmentation;
  DomainSegmentArray<T, D> domainSegments;

  // std::vector<SizeType> activePointIdOffsets;
  std::vector<SizeType> pointIdOffsets;

  // Friend classes
  // Iterators need private access to get the values
  template <class> friend class BaseIterator;
  template <class> friend class SparseIterator;
  template <class> friend class SparseOffsetIterator;

  template <class S, class Type> static S &writeValue(S &stream, Type &value) {
    if constexpr (std::is_arithmetic_v<Type>) {
      stream.write(reinterpret_cast<const char *>(&value), sizeof(value));
    } else {
      value->serialize(stream);
    }
    return stream;
  }

  template <class S, class Type> static S &readValue(S &stream, Type &value) {
    if constexpr (std::is_arithmetic_v<Type>) {
      stream.read(reinterpret_cast<char *>(&value), sizeof(value));
    } else {
      value->deserialize(stream);
    }
    return stream;
  }

public:
  // CONSTRUCTORS
  Domain() = default;
  Domain(const Domain &) = delete;

  // create empty level set with one undefined run
  explicit Domain(Grid<D> &g, SizeType runType = RunTypeValues::UNDEF_PT)
      : grid(&g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
  };

  explicit Domain(Grid<D> *g, SizeType runType = RunTypeValues::UNDEF_PT)
      : grid(g) {
    initialize();
    domainSegments[0].insertNextUndefinedRunType(grid->getMinIndex(), runType);
  };

  // create empty level set with one undefined value
  Domain(Grid<D> &g, T value) : grid(&g) {
    initialize();
    domainSegments[0].insertNextUndefinedPoint(grid->getMinIndex(), value);
  };

  Domain(Grid<D> *g, T value) : grid(g) {
    initialize();
    domainSegments[0].insertNextUndefinedPoint(grid->getMinIndex(), value);
  };

  void deepCopy(const Domain &passedDomain) {
    assert(grid == &(passedDomain.getGrid()));
    pointIdOffsets = passedDomain.pointIdOffsets;
    segmentation = passedDomain.segmentation;
    domainSegments.deepCopy(grid, passedDomain.domainSegments);
  }

  void deepCopy(Grid<D> *passedGrid, const Domain &passedDomain) {
    deepCopy(*passedGrid, passedDomain);
  }

  void deepCopy(Grid<D> &passedGrid, const Domain &passedDomain) {
    grid = &passedGrid;
    pointIdOffsets = passedDomain.pointIdOffsets;
    segmentation = passedDomain.segmentation;
    domainSegments.deepCopy(grid, passedDomain.domainSegments);
  }

  void shallowCopy(const Domain &passedDomain) {
    grid = passedDomain.grid;
    pointIdOffsets = passedDomain.pointIdOffsets;
    domainSegments.shallowCopy(passedDomain.domainSegments);
    segmentation = passedDomain.segmentation;
  }

  // Inline member functions
  void print(std::ostream &out = std::cout) const {
    for (SizeType i = 0; i != domainSegments.getNumberOfSegments(); ++i) {
      domainSegments[i].print(out);
      out << std::endl;
    }
  }

  const Grid<D> &getGrid() const { return *grid; }

  Grid<D> &getGrid() { return *grid; }

  DomainSegmentType &getDomainSegment(unsigned i) { return domainSegments[i]; }

  const DomainSegmentType &getDomainSegment(unsigned i) const {
    return domainSegments[i];
  }

  SizeType getPointIdOffset(unsigned i) const { return pointIdOffsets[i]; }

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
    return static_cast<unsigned>(domainSegments.getNumberOfSegments());
  }

  unsigned getNumberOfRuns(int segmentId, int dimension) const {
    return domainSegments[segmentId].getNumberOfRuns(dimension);
  }

  const IndexPoints &getSegmentation() const { return segmentation; }

  IndexType getMinRunBreak(int dim) const {
    IndexType minBreak = domainSegments[0].getMinRunBreak(dim);
    for (unsigned i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      minBreak = std::min(minBreak, domainSegments[i].getMinRunBreak(dim));
    }
    return minBreak;
  }

  IndexType getMaxRunBreak(int dim) const {
    IndexType maxBreak = domainSegments[0].getMaxRunBreak(dim);
    for (unsigned i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      maxBreak = std::max(maxBreak, domainSegments[i].getMaxRunBreak(dim));
    }
    return maxBreak;
  }

  Index<D> getMinRunBreak() const {
    Index<D> minBreak;
    for (unsigned i = 0; i < D; ++i) {
      minBreak[i] = getMinRunBreak(i);
    }
    return minBreak;
  }

  Index<D> getMaxRunBreak() const {
    Index<D> maxBreak;
    for (unsigned i = 0; i < D; ++i) {
      maxBreak[i] = getMaxRunBreak(i);
    }
    return maxBreak;
  }

  void getDomainBounds(IndexType *bounds) {
    for (unsigned i = 0; i < D; ++i) {
      if (grid->getBoundaryConditions(i) == BoundaryType::INFINITE_BOUNDARY ||
          grid->getBoundaryConditions(i) ==
              BoundaryType::NEG_INFINITE_BOUNDARY) {
        bounds[2 * i] = getMinRunBreak(i);
      } else {
        bounds[2 * i] = grid->getMinGridPoint(i);
      }

      if (grid->getBoundaryConditions(i) == BoundaryType::INFINITE_BOUNDARY ||
          grid->getBoundaryConditions(i) ==
              BoundaryType::POS_INFINITE_BOUNDARY) {
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
  void insertNextUndefinedRunType(int sub, const V &point, SizeType rt) {
    // this function sets the sign of an undefined run starting at indices
    // "point" the indices of the new grid point "indices" must be greater
    // (lexiographically greater) than the indices of the current "last" grid
    // point
    domainSegments[sub].insertNextUndefinedRunType(point, rt);
  }

  template <class V>
  void insertNextUndefinedPoint(int sub, const V &point, ValueType value) {
    // this function sets the sign of an undefined run starting at indices
    // "point" the indices of the new grid point "indices" must be greater
    // (lexiographically greater) than the indices of the current "last" grid
    // point
    domainSegments[sub].insertNextUndefinedPoint(point, value);
  }

  // create new segmentation and domainSegments
  void initialize(
      const IndexPoints &p = IndexPoints(),
      const AllocationType<SizeType, D> &a = AllocationType<SizeType, D>()) {

    segmentation = p;

    domainSegments.initialize(segmentation.size() + 1, *grid,
                              a); // clear old segments and create new ones

    for (SizeType i = 1; i < domainSegments.getNumberOfSegments(); ++i) {
      DomainSegment<T, D> &s = domainSegments[i];

      s.insertNextUndefinedRunType(grid->getMinIndex(),
                                   grid->decrementIndices(segmentation[0]),
                                   RunTypeValues::SEGMENT_PT);

      for (SizeType j = 1; j < i; ++j) {
        s.insertNextUndefinedRunType(segmentation[j - 1],
                                     grid->decrementIndices(segmentation[j]),
                                     RunTypeValues::SEGMENT_PT + j);
      }
    }
  }

  // distribute data points evenly across DomainSegments and add SEGMENT_PT
  // as boundary markers
  void finalize() {

    for (SizeType i = 0; i + 1 < domainSegments.getNumberOfSegments(); ++i) {

      DomainSegment<T, D> &s = domainSegments[i];

      for (SizeType j = i + 1; j < segmentation.size(); ++j) {
        s.insertNextUndefinedRunType(segmentation[j - 1],
                                     grid->decrementIndices(segmentation[j]),
                                     RunTypeValues::SEGMENT_PT + j);
      }
      s.insertNextUndefinedRunType(
          segmentation.back(), grid->getMaxGridPoint(),
          RunTypeValues::SEGMENT_PT +
              static_cast<SizeType>(segmentation.size()));
    }

    // TODO: for now do not save pointId offsets as they can easily be found
    // from cumulating domainSegments[i].getNumberOfPoints()
    // calculate id-offsets
    pointIdOffsets.clear();
    // active_pointIdOffsets.clear();
    pointIdOffsets.push_back(0);
    // activePointId_offsets.push_back(0);
    for (SizeType i = 0; i < domainSegments.getNumberOfSegments() - 1; ++i) {
      pointIdOffsets.push_back(pointIdOffsets.back() +
                               domainSegments[i].getNumberOfPoints());
      //     activePointId_offsets.push_back(activePointId_offsets.back()+domainSegments[i].num_active_pts());
    }

    // assert(activePointId_offsets.size()==DomainSegment.size());
    assert(pointIdOffsets.size() == domainSegments.getNumberOfSegments());
  }

  /// Converts a pointId(given by lexicographical order of points) to a spatial
  /// coordinate
  Index<D> ptIdToCoordinate(SizeType pt) const {
    Index<D> pointCoords;

    // find domainSegment of point
    const int segment = int(
        std::upper_bound(pointIdOffsets.begin() + 1, pointIdOffsets.end(), pt) -
        (pointIdOffsets.begin() + 1));

    pt -= pointIdOffsets[segment]; // local point id

    const DomainSegmentType &s = domainSegments[segment];

    // find right PointID by bisection
    for (int g = 0; g < D; ++g) {
      SizeType min = 0;
      SizeType max = static_cast<SizeType>(s.runTypes[g].size()) - 1;

      while (!s.isPtIdDefined(s.runTypes[g][min]))
        ++min;
      while (!s.isPtIdDefined(s.runTypes[g][max]))
        --max;

      while (min != max) {
        SizeType mid = (max + min + 1) / 2;
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

      pt =
          static_cast<SizeType>(std::upper_bound(s.startIndices[g].begin() + 1,
                                                 s.startIndices[g].end(), min) -
                                (s.startIndices[g].begin() + 1));

      pointCoords[g] += s.getRunStartCoord(g, pt, min);
    }

    return pointCoords;
  }

  /// finds the ideal coordinates at which to break between domainSegments
  /// for balanced distribution of points
  IndexPoints getNewSegmentation() const {
    IndexPoints tempSegmentation;

    int n = 1;
#ifdef _OPENMP
    n = omp_get_max_threads();
#endif
    SizeType n_pts = getNumberOfPoints(); // number of defined grid points
    SizeType sum = 0;

    for (unsigned int g = 0; g < static_cast<unsigned int>(n - 1); ++g) {
      sum += n_pts / n + ((n_pts % n) > g);
      // TODO: is that if really necessary?
      if (sum != n_pts)
        tempSegmentation.push_back(ptIdToCoordinate(sum));
    }

    return tempSegmentation;
  }

  /// allocation_type allocates the required sizes to num_values and num_runs
  /// for all the domainSegments members num_values[0] is to contain level set
  /// values, num_values[i] contains the start indices at the i-th dimension
  /// num_runs[i] is to contain the run types at the i-th dimension
  AllocationType<SizeType, D> getAllocation() const {
    AllocationType<SizeType, D> a;
    for (unsigned i = 0; i < domainSegments.getNumberOfSegments(); ++i) {
      AllocationType<SizeType, D> b = domainSegments[i].getAllocation();
      a.num_values = Max(a.num_values, b.num_values);
      a.num_runs = Max(a.num_runs, b.num_runs);
    }

    return a * domainSegments.getNumberOfSegments();
  }

  /// distribute points evenly across domainSegments, so that they can be
  /// iterated over by separate threads efficiently
  void segment() {
    if (getNumberOfPoints() == 0) {
      return;
    }
    Domain newDomain(grid);
    newDomain.initialize(getNewSegmentation(), getAllocation());

#pragma omp parallel num_threads(                                              \
        int(newDomain.domainSegments.getNumberOfSegments()))
    {
      int p = 0;
#ifdef _OPENMP
      p = omp_get_thread_num();
#endif

      DomainSegmentType &s = newDomain.domainSegments[p];

      Index<D> startOfSegment =
          (p == 0) ? grid->getMinIndex() : newDomain.segmentation[p - 1];

      Index<D> endOfSegment =
          (p != static_cast<int>(newDomain.segmentation.size()))
              ? newDomain.segmentation[p]
              : grid->getMaxGridPoint();

      for (SparseIterator<Domain> it(*this, startOfSegment);
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
      Domain newDomain(grid);
      newDomain.initialize(); // initialize with only one segment

      DomainSegmentType &s = newDomain.domainSegments[0];
      for (ConstSparseIterator<Domain> it(*this); !it.isFinished(); ++it) {
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
    using namespace hrleUtil;
    // write identifier
    stream << "hrleDomain";

    bool structureIsSerial = true;
    if (getNumberOfSegments() > 1) {
      desegment();
      structureIsSerial = false;
    }

    // get domain bounds
    IndexType bounds[2 * D];
    getDomainBounds(bounds);

    // START INDICES, RUN TYPES, RUN BREAKS FOR EACH DIMENSION
    for (int dim = D - 1; dim >= 0; --dim) {
      // get start indices, runbreaks and runtypes
      const std::vector<SizeType> &startIndices =
          domainSegments[0].startIndices[dim];
      const std::vector<SizeType> &runTypes = domainSegments[0].runTypes[dim];
      const std::vector<IndexType> &runBreaks =
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
        auto numberOfValues = uint32_t(startIndices.size());
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
      std::vector<SizeType>
          definedRunIndices; // store all indices of defined runtypes
      if (bitsPerRunType > 0) {
        int count = 8 / bitsPerRunType - 1;
        unsigned char byte = 0;

        // each runType needs at least one byte
        if (bitsPerRunType > 4) {
          for (unsigned long runType : runTypes) {
            SizeType PtId = 0;
            // if undefined point, need to shift id
            if (runType >= RunTypeValues::UNDEF_PT)
              PtId = runType - RunTypeValues::UNDEF_PT + 1;
            else
              definedRunIndices.push_back(runType);

            stream.write((char *)&PtId, (bitsPerRunType - 1) / 8 + 1);
          }
        } else { // can fit more than one value in a byte
          for (unsigned long runType : runTypes) {
            SizeType PtId = 0;
            if (runType >= RunTypeValues::UNDEF_PT)
              PtId = runType - RunTypeValues::UNDEF_PT + 1;
            else
              definedRunIndices.push_back(runType);

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
      }

      // Write indices of defined runtypes; only save the difference to the next
      // defined runtype write the first runtype(always 0) explicitly; makes
      // reading easier
      if (!definedRunIndices.empty()) {
        stream.write((char *)&definedRunIndices[0], bytesPerIndex);
        for (unsigned int i = 0; i < definedRunIndices.size() - 1; i++) {
          unsigned long diff = definedRunIndices[i + 1] - definedRunIndices[i];
          stream.write((char *)&diff, bytesPerIndex);
        }
      }

      // Write runbreaks
      for (int runBreak : runBreaks) {
        stream.write((char *)&runBreak, bytesPerRunBreak);
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
      std::vector<ValueType> &definedValues = domainSegments[0].definedValues;
      for (typename std::vector<ValueType>::const_iterator it =
               definedValues.begin();
           it != definedValues.end(); ++it) {
        writeValue(stream, *it);
        // stream.write(reinterpret_cast<char*>*it) //<< *it;
      }

      std::vector<ValueType> &undefinedValues =
          domainSegments[0].undefinedValues;
      for (typename std::vector<ValueType>::const_iterator it =
               undefinedValues.begin();
           it != undefinedValues.end(); ++it) {
        writeValue(stream, *it);
        // stream << *it;
      }
    }

    // if the Domain was segmented before, segment it again
    if (!structureIsSerial) {
      segment();
    }

    return stream;
  }

  std::istream &deserialize(std::istream &stream) {
    // check identifier
    char identifier[11] = {}; // 1 more for string constructor
    stream.read(identifier, 10);
    if (std::string(identifier).compare(0, 10, "hrleDomain")) {
      std::cout
          << "Reading Domain from stream failed. Header could not be found."
          << std::endl;
      return stream;
    }

    // READ HRLE PROPERTIES
    for (int dim = D - 1; dim >= 0; --dim) {
      // get the start indices, runtypes and runbreaks vectors
      std::vector<SizeType> &startIndices = domainSegments[0].startIndices[dim];
      std::vector<SizeType> &runTypes = domainSegments[0].runTypes[dim];
      std::vector<IndexType> &runBreaks = domainSegments[0].runBreaks[dim];

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
        if (numberOfStartIndices > 0) {
          startIndices.push_back(0);
          for (unsigned i = 1; i < numberOfStartIndices; ++i) {
            unsigned long long current = 0;
            stream.read((char *)&current, bytesPerIndex);
            sum += current;
            startIndices.push_back(SizeType(sum));
          }
        }
      }

      // READ RUN TYPES
      if (bitsPerRunType > 0) {
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
            SizeType current = 0;
            stream.read((char *)&current, (bitsPerRunType - 1) / 8 + 1);
            if (current == 0) { // defined run
              // for a defined run, store 0, so it can be replaced later
              runTypes[i] = 0;
            } else { // undefined run
              runTypes[i] = current - 1 + RunTypeValues::UNDEF_PT;
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
              SizeType current = (byte & bitMask) >> (8 - bitsPerRunType);
              byte <<= bitsPerRunType;
              if (current == 0) { // defined run
                runTypes[i * numberOfValuesPerByte + j] = 0;
              } else { // undefined run
                runTypes[i * numberOfValuesPerByte + j] =
                    current - 1 + RunTypeValues::UNDEF_PT;
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
          runBreaks.push_back(static_cast<IndexType>(runBreak));
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
        ValueType definedValue;
        readValue(stream, definedValue);
        domainSegments[0].definedValues.push_back(definedValue);
      }

      domainSegments[0].undefinedValues.clear();
      // UNDEFINED VALUES
      for (unsigned i = 0; i < numberOfUndefinedValues; ++i) {
        ValueType undefinedValue;
        readValue(stream, undefinedValue);
        domainSegments[0].undefinedValues.push_back(undefinedValue);
      }
    }

    finalize();
    segment();

    return stream;
  }
};
} // namespace viennahrle

#endif // HRLE_DOMAIN_HPP_
