#ifndef HRLE_DOMAIN_SEGMENT_HPP
#define HRLE_DOMAIN_SEGMENT_HPP

#include <cassert>
#include <iomanip>
#include <sstream>

#include "hrleAllocationType.hpp"
#include "hrleGrid.hpp"
#include "hrleRunTypeValues.hpp"
#include "hrleTypes.hpp"

namespace viennahrle {
template <class T = double, int D = 3> class DomainSegment {

  const Grid<D> *grid;

public:
  // TYPEDEFS
  typedef T ValueType;

  // the following vectors are used to store the run-length-encoded data
  // structure for more details on that data structure see B. Houston, M.B.
  // Nielsen, C. Batty, O. Nilsson, K. Museth, "Hierarchical RLE-Level Set -
  // A Compact and Versatile Deformable Surface Representation", ACM Trans.
  // Graph. 25/1, pp. 151-175, 2006.
  std::vector<SizeType> startIndices[D];
  std::vector<SizeType> runTypes[D];
  std::vector<IndexType> runBreaks[D];
  std::vector<ValueType>
      definedValues; // this vector keeps the values of all defined
                     // grid points its size is therefore equal to the number of
                     // defined grid points

  std::vector<ValueType> undefinedValues; // this vector keeps all the values of
                                          // undefined runs ( background values)

  // std::vector<SizeType>
  //     activePointIds; // an additional vector, which has the same size as the
  // distances-vector in case the defined grid point is
  // active in terms of the sparse field level set method
  // (the level set value is in the range [-0.5,0.5]) this
  // vector keeps the active point id the active point ids
  // are given to the active grid points in lexicographical
  // order starting with index 0 therefore the largest
  // active point ID equals the number of active grid
  // points-1 if the grid point is not active, its value is
  // set to the constant RunTypeValues::INACTIVE_PT

  // SizeType numberOfActivePoints; // numberOfActivePoints stores the number
  // of active points

  // STATIC FUNCTIONS
  static bool isPtIdDefined(const SizeType r) { // returns if the grid point
                                                // given by the "getPointId"
    return (r < RunTypeValues::UNDEF_PT);       // is a defined grid point
  }

  // DYNAMIC FUNCTIONS
  AllocationType<SizeType, D> getAllocation() const {
    // AllocationType allocates the required sizes to num_values and
    // num_runs num_values[0] is to contain defined values, num_values[i]
    // contains the start indices at the i-th dimension num_runs[i] is to
    // contain the run types at the i-th dimension
    AllocationType<SizeType, D> a;
    a.num_values[0] = static_cast<SizeType>(definedValues.size());
    a.num_runs[0] = static_cast<SizeType>(runTypes[0].size());
    for (int i = 1; i < D; ++i) {
      a.num_values[i] = static_cast<SizeType>(startIndices[i - 1].size());
      a.num_runs[i] = static_cast<SizeType>(runTypes[i].size());
    }

    // Check that the vector sizes make sense for the H-RLE structure
    assert(runBreaks[D - 1].size() == a.num_runs[D - 1] - 1);
    for (int i = 1; i < D; ++i) {
      assert(runBreaks[i - 1].size() == a.num_runs[i - 1] - a.num_values[i]);
    }
    return a;
  }

  // No default ctor
  DomainSegment() = delete;

  /// Constructor for the DomainSegment class
  /// Takes a pointer to a grid_type and an AllocationType to reserve the
  /// memory required for the start indexes, run types and run breaks arrays, as
  /// well as the LS values
  DomainSegment(const Grid<D> &g, const AllocationType<SizeType, D> &a)
      : grid(&g) {
    definedValues.reserve(a.num_values[0]);
    runTypes[0].reserve(a.num_runs[0]);

    for (int i = 1; i < D; ++i) {
      startIndices[i - 1].reserve(a.num_values[i]);
      runBreaks[i - 1].reserve(a.num_runs[i - 1] - a.num_values[i]);
      runTypes[i].reserve(a.num_runs[i]);
    }

    startIndices[D - 1].push_back(0);
    runBreaks[D - 1].reserve(a.num_runs[D - 1] - 1);
  };

  DomainSegment(const Grid<D> &g, const DomainSegment &s) : grid(&g) {
    for (int i = 0; i < D; ++i)
      startIndices[i] = s.startIndices[i];
    for (int i = 0; i < D; ++i)
      runTypes[i] = s.runTypes[i];
    for (int i = 0; i < D; ++i)
      runBreaks[i] = s.runBreaks[i];
    definedValues = s.definedValues;
    undefinedValues = s.undefinedValues;
    grid = &g;
  };

  DomainSegment &operator=(const DomainSegment &s) {
    if (this == &s) // handle self-assignment
      return *this;

    for (int i = 0; i < D; ++i)
      startIndices[i] = s.startIndices[i];
    for (int i = 0; i < D; ++i)
      runTypes[i] = s.runTypes[i];
    for (int i = 0; i < D; ++i)
      runBreaks[i] = s.runBreaks[i];
    definedValues = s.definedValues;
    undefinedValues = s.undefinedValues;
    // activePointIds = s.activePointIds;
    // numberOfActivePoints = s.numberOfActivePoints;
    grid = s.grid;
    return *this;
  }

  /// returns the start index of the run given by startIndicesPos and
  /// runTypePos see the definition of the HRLE-data structure for more
  /// details
  IndexType getRunStartCoord(int dim, SizeType startIndicesPos,
                             SizeType runTypePos) const {

    if (runTypePos == startIndices[dim][startIndicesPos]) {
      return grid->getMinIndex(dim);
    } else {
      return runBreaks[dim][runTypePos - startIndicesPos - 1];
    }
  }

  /// returns the end index of the run given by startIndicesPos and
  /// runTypePos NOTE: the end index is not included by the run see the
  /// definition of the HRLE-data structure for more details
  IndexType getRunEndCoord(int dim, SizeType startIndicesPos,
                           SizeType runTypePos) const {
    if (runTypePos + 1 < getStartIndex(dim, startIndicesPos + 1)) {
      return runBreaks[dim][runTypePos - startIndicesPos] - 1;
    } else {
      return grid->getMaxGridPoint(dim);
    }
  }

  /// returns the starting index of the runTypes array
  /// if startIndicesPos is equal to the size of startIndices array
  /// the size of the run_types array is returned
  SizeType getStartIndex(int dim, SizeType startIndicesPos) const {
    if (startIndicesPos == startIndices[dim].size())
      return static_cast<SizeType>(runTypes[dim].size());
    return startIndices[dim][startIndicesPos];
  }

  IndexType getRunBreak(int dim,
                        int runbreak = std::numeric_limits<int>::max()) const {
    if (runBreaks[dim].size() == 0)
      return 0;
    if (runbreak == std::numeric_limits<int>::max())
      return runBreaks[dim].back();
    return runBreaks[dim][runbreak];
  }

  IndexType getMaxRunBreak(int dim) const {
    if (runBreaks[dim].size() == 0)
      return 0;
    return *std::max_element(runBreaks[dim].begin(), runBreaks[dim].end());
  }

  IndexType getMinRunBreak(int dim) const {
    if (runBreaks[dim].size() == 0)
      return 0;
    return *std::min_element(runBreaks[dim].begin(), runBreaks[dim].end());
  }

  /// this function gives an estimation of the used memory of the level
  /// set function. however, the allocated memory is much higher because the
  /// STL-vector allocates more memory than needed, to accelerate
  /// push_back-operations the allocated memory can be obtained by the
  /// "allocated_memory"-function defined below NOTE: since this function is
  /// based on  the sizeof-operator it cannot take the memory needed for the
  ///      GridTraitsType (see grid.hpp) into account if dynamic data
  ///      structures are used for example Furthermore due to data structure
  ///      alignment the result can be somehow inaccurate
  unsigned long int getMemoryFootprint() const {
    unsigned long int x = sizeof(DomainSegment);

    for (int i = 0; i < D; i++) {
      x += sizeof(SizeType) * startIndices[i].size();
      x += sizeof(SizeType) * runTypes[i].size();
      x += sizeof(IndexType) * runBreaks[i].size();
    }
    x += sizeof(ValueType) * definedValues.size();
    // x += sizeof(SizeType) * activePointIds.size();

    return x;
  }

  /// this function gives an estimation of the allocated memory for the level
  /// set function.
  /// NOTE: since this function is based on  the sizeof-operator it cannot take
  /// the memory needed for the
  ///      GridTraitsType (see grid.hpp) into account if dynamic data
  ///      structures are used for example Furthermore due to data structure
  ///      alignment the result can be somehow inaccurate
  unsigned long int getAllocatedMemory() const {
    unsigned long int x = sizeof(DomainSegment);

    for (int i = 0; i < D; i++) {
      x += sizeof(SizeType) * startIndices[i].capacity();
      x += sizeof(SizeType) * runTypes[i].capacity();
      x += sizeof(IndexType) * runBreaks[i].capacity();
    }
    x += sizeof(ValueType) * definedValues.capacity();
    // x += sizeof(SizeType) * activePointIds.capacity();

    return x;
  }

  /// this function returns the number of defined grid points
  SizeType getNumberOfPoints() const {
    return static_cast<SizeType>(definedValues.size());
  };

  /// this function returns the number of distinct undefined values
  SizeType getNumberOfUndefinedValues() const {
    return static_cast<SizeType>(undefinedValues.size());
  }

  /// returns the number of distinct runs in the dimension <dimension>
  /// for negative <dimension> it returns the number of defined points
  SizeType getNumberOfRuns(int dimension) const {
    if (dimension >= 0)
      return runTypes[dimension].size();
    return definedValues.size();
  }

  template <class V>
  void insertNextUndefinedRunType(V start_point, const V &end_point,
                                  SizeType rt) {

    if (start_point > end_point)
      return; // in this case, do not add the point

    for (int dim = 0; dim < D - 1; ++dim) {
      if (start_point[dim] != grid->getMinGridPoint(dim)) {
        if (start_point[dim] <= grid->getMaxGridPoint(dim)) {
          insertNextUndefinedRunType(start_point, rt);
        }

        start_point[dim] = grid->getMinIndex(dim);
        ++start_point[dim + 1];
      }

      if (start_point > end_point)
        return;
    }

    if (start_point[D - 1] <= grid->getMaxGridPoint(D - 1)) {
      insertNextUndefinedRunType(start_point, rt);
    }
  }

  template <class V>
  void insertNextUndefinedRunType(const V &point, SizeType rt) {

    int level;
    for (level = 0; level < D; ++level) {
      if (point[level] != grid->getMinIndex(level))
        break;
    }

    SizeType old_sign = 0;
    int dim;

    for (dim = D - 1; dim > level; --dim) {

      if (runTypes[dim].size() ==
          startIndices[dim].back()) { // if there is no run
        if (point[dim] != grid->getMinIndex(dim)) {
          runTypes[dim].push_back(old_sign);
          runBreaks[dim].push_back(point[dim]);
        }
        runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
        startIndices[dim - 1].push_back(SizeType(runTypes[dim - 1].size()));
      } else if (!isPtIdDefined(
                     runTypes[dim].back())) { // if there is a defined run
        old_sign = runTypes[dim].back();
        if (old_sign == rt)
          return;
        if (runTypes[dim].size() ==
            startIndices[dim].back() + 1) { // if there is a single run
          if (point[dim] == grid->getMinIndex(dim)) {
            runTypes[dim].back() = SizeType(startIndices[dim - 1].size());
          } else {
            runBreaks[dim].push_back(point[dim]);
            runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
          }
        } else { // if there are more than one runs
          if (point[dim] == runBreaks[dim].back()) {
            runTypes[dim].pop_back();
            if (!isPtIdDefined(runTypes[dim].back())) {
              runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
            } else
              runBreaks[dim].pop_back();
          } else {
            runBreaks[dim].push_back(point[dim]);
            runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
          }
        }
        startIndices[dim - 1].push_back(SizeType(runTypes[dim - 1].size()));
      }
    }

    if (runTypes[dim].size() ==
        startIndices[dim].back()) { // if there is no run
      if (point[dim] != grid->getMinIndex(dim)) {
        runTypes[dim].push_back(old_sign);
        runBreaks[dim].push_back(point[dim]);
      }
      runTypes[dim].push_back(rt);
    } else if (!isPtIdDefined(
                   runTypes[dim].back())) { // if there is an defined run
      old_sign = runTypes[dim].back();
      if (old_sign == rt)
        return;
      if (runTypes[dim].size() ==
          startIndices[dim].back() + 1) { // if there is a single run
        if (point[dim] == grid->getMinIndex(dim)) {
          runTypes[dim].back() = rt;
        } else {
          runBreaks[dim].push_back(point[dim]);
          runTypes[dim].push_back(rt);
        }
      } else { // if there are more than one runs
        if (point[dim] == runBreaks[dim].back()) {
          runTypes[dim].back() = rt;
        } else {
          runBreaks[dim].push_back(point[dim]);
          runTypes[dim].push_back(rt);
        }
      }
    } else {
      runBreaks[dim].push_back(point[dim]);
      runTypes[dim].push_back(rt);
    }
  }

  template <class V>
  void insertNextUndefinedPoint(const V &startPoint, const V &endPoint,
                                ValueType value) {
    // if undefined value already exists, use its runtype,
    // if it does not use the next available undefined runtype
    auto it = std::find(undefinedValues.begin(), undefinedValues.end(), value);
    SizeType runType;
    if (it != undefinedValues.end()) {
      runType =
          RunTypeValues::UNDEF_PT +
          std::distance(
              undefinedValues.begin(),
              it); // undefined runtypes start at RunTypeValues::UNDEF_PT+1
    } else {
      runType = RunTypeValues::UNDEF_PT + undefinedValues.size();
      undefinedValues.push_back(value);
    }

    insertNextUndefinedRunType(startPoint, endPoint, runType);
  }

  template <class V>
  void insertNextUndefinedPoint(const V &point, ValueType value) {
    // if undefined value already exists, use its runtype,
    // if it does not, use the next available undefined runtype
    auto it = std::find(undefinedValues.begin(), undefinedValues.end(), value);
    SizeType runType;
    if (it != undefinedValues.end()) {
      runType = RunTypeValues::UNDEF_PT +
                SizeType(std::distance(undefinedValues.begin(), it));
    } else {
      runType = RunTypeValues::UNDEF_PT + SizeType(undefinedValues.size());
      undefinedValues.push_back(value);
    }

    insertNextUndefinedRunType(point, runType);
  }

  /// Inserts an undefined point into the HRLE structure.
  /// CAREFUL: If the same point is inserted twice, the structure might break!
  template <class V>
  void insertNextDefinedPoint(const V &point, ValueType value) {

    int level;
    for (level = 0; level < D; ++level) {
      if (point[level] != grid->getMinIndex(level))
        break;
    }

    SizeType old_sign = 0;

    for (int dim = D - 1; dim > 0; --dim) {
      if (runTypes[dim].size() ==
          startIndices[dim].back()) { // if there is no run
        if (point[dim] != grid->getMinIndex(dim)) {
          runTypes[dim].push_back(old_sign);
          runBreaks[dim].push_back(point[dim]);
        }
        runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
        startIndices[dim - 1].push_back(SizeType(runTypes[dim - 1].size()));
      } else if (!isPtIdDefined(
                     runTypes[dim].back())) { // if there is an undefined run
        old_sign = runTypes[dim].back();
        if (runTypes[dim].size() ==
            startIndices[dim].back() + 1) { // if there is a single run
          if (point[dim] == grid->getMinIndex(dim)) {
            runTypes[dim].back() = SizeType(startIndices[dim - 1].size());
          } else {
            runBreaks[dim].push_back(point[dim]);
            runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
          }
        } else { // if there are more than one runs
          if (point[dim] == runBreaks[dim].back()) {
            runTypes[dim].pop_back();
            if (!isPtIdDefined(runTypes[dim].back())) {
              runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
            } else {
              runBreaks[dim].pop_back();
            }
          } else {
            runBreaks[dim].push_back(point[dim]);
            runTypes[dim].push_back(SizeType(startIndices[dim - 1].size()));
          }
        }
        startIndices[dim - 1].push_back(SizeType(runTypes[dim - 1].size()));
      } else { // it is a defined run
        if (dim <= level)
          startIndices[dim - 1].push_back(SizeType(runTypes[dim - 1].size()));
      }
    }

    if (runTypes[0].size() == startIndices[0].back()) { // if there is no run
      if (point[0] != grid->getMinIndex(0)) {
        runTypes[0].push_back(old_sign);
        runBreaks[0].push_back(point[0]);
      }
      runTypes[0].push_back(SizeType(definedValues.size()));
    } else if (!isPtIdDefined(
                   runTypes[0].back())) { // if there is an undefined run
      old_sign = runTypes[0].back();
      if (runTypes[0].size() ==
          startIndices[0].back() + 1) { // if there is a single run
        if (point[0] == grid->getMinIndex(0)) {
          runTypes[0].back() = SizeType(definedValues.size());
        } else {
          runBreaks[0].push_back(point[0]);
          runTypes[0].push_back(SizeType(definedValues.size()));
        }
      } else { // if there are more than one runs
        if (point[0] == runBreaks[0].back()) {
          runTypes[0].pop_back();
          if (!isPtIdDefined(runTypes[0].back())) {
            runTypes[0].push_back(SizeType(definedValues.size()));
          } else {
            runBreaks[0].pop_back();
          }
        } else {
          runBreaks[0].push_back(point[0]);
          runTypes[0].push_back(SizeType(definedValues.size()));
        }
      }
    }

    definedValues.push_back(value);

    // if (abs(value) <= ValueType(0.5)) {
    //   activePointIds.push_back(numberOfActivePoints);
    //   ++numberOfActivePoints;
    //
    // } else {
    //   activePointIds.push_back(RunTypeValues::INACTIVE_PT);
    // }
  }

  void print(std::ostream &out = std::cout) const {
    std::ostringstream oss;

    out << std::endl;
    out << std::string(20, '-') << " HRLE Data Structure "
        << std::string(20, '-') << std::endl
        << std::endl;

    for (int dim = D - 1; dim >= 0; --dim) {
      out << dim << " startIndices: " << startIndices[dim].size();
      int c = 0;
      for (auto it = startIndices[dim].begin(); it != startIndices[dim].end();
           ++it) {
        if (c % 10 == 0)
          out << std::endl;
        ++c;
        out << std::setw(8) << *it;
      }
      out << std::endl;

      out << dim << " run_types: " << runTypes[dim].size();
      c = 0;
      for (auto it = runTypes[dim].begin(); it != runTypes[dim].end(); ++it) {
        if (c % 10 == 0)
          out << std::endl;
        ++c;

        if ((*it) >= RunTypeValues::SEGMENT_PT) {
          out << std::setw(8) << "SEG" << (*it) - RunTypeValues::SEGMENT_PT;
        } else if ((*it) >= RunTypeValues::UNDEF_PT) {
          out << std::setw(8) << "U" << (*it) - RunTypeValues::UNDEF_PT;
        } else {
          out << std::setw(8) << (*it);
        }
      }
      out << std::endl;

      out << dim << " run_breaks: " << runBreaks[dim].size();
      c = 0;
      for (auto it = runBreaks[dim].begin(); it != runBreaks[dim].end(); ++it) {
        if (c % 10 == 0)
          out << std::endl;
        ++c;
        out << std::setw(8) << *it;
      }
      out << std::endl;
    }

    out << "definedValues: " << definedValues.size();
    int c = 0;
    for (auto it = definedValues.begin(); it != definedValues.end(); ++it) {
      if (c % 10 == 0)
        out << std::endl;
      ++c;
      out << std::setw(8) << std::fixed << *it << " ";
    }
    out << std::endl;
    out << std::endl;

    out << "undefinedValues: " << undefinedValues.size();
    for (unsigned i = 0; i < undefinedValues.size(); ++i) {
      if (i % 10 == 0)
        out << std::endl;
      out << std::setw(16) << std::defaultfloat << undefinedValues[i];
    }
    out << std::endl;
    out << std::endl;

    out << std::string(60, '-') << std::endl << std::endl;
  }
};
} // namespace viennahrle

#endif // HRLE_DOMAIN_SEGMENT_HPP
