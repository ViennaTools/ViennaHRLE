#ifndef HRLE_RUN_TYPE_VALUES_HPP
#define HRLE_RUN_TYPE_VALUES_HPP

#include <limits>

#include "hrleTypes.hpp"

// the number of possible undefined run value, ergo how many different
// background values there can be
#ifndef MAX_NUMBER_OF_BACKGROUND_VALUES
#define MAX_NUMBER_OF_BACKGROUND_VALUES 10000
#endif // MAX_NUMBER_OF_BACKGROUND_VALUES

// the maximum number of threads allowed
#ifndef MAX_NUMBER_OF_OMP_THREADS
#define MAX_NUMBER_OF_OMP_THREADS 100
#endif // MAX_NUMBER_OF_OMP_THREADS

namespace viennahrle {
struct RunTypeValues {
  static constexpr SizeType SEGMENT_PT =
      std::numeric_limits<SizeType>::max() -
      (MAX_NUMBER_OF_OMP_THREADS); // this value signifies end the of
                                   // the current hrleDomainSegment

  static constexpr SizeType UNDEF_PT = // undef runs have ids up to SEGMENT_PT
      std::numeric_limits<SizeType>::max() -
      (MAX_NUMBER_OF_BACKGROUND_VALUES + MAX_NUMBER_OF_OMP_THREADS);

  static constexpr SizeType INACTIVE_PT = std::numeric_limits<SizeType>::max();
};
} // namespace viennahrle

#endif // HRLE_RUN_TYPE_VALUES_HPP
