#ifndef HRLE_RUN_TYPE_VALUES_HPP
#define HRLE_RUN_TYPE_VALUES_HPP

#include <limits>

#include "hrleSizeType.hpp"

// the number of possible undefined run value, ergo how many different
// background values there can be
#ifndef MAX_NUMBER_OF_BACKGROUND_VALUES
#define MAX_NUMBER_OF_BACKGROUND_VALUES 10000
#endif // MAX_NUMBER_OF_BACKGROUND_VALUES

// the maximum number of threads allowed
#ifndef MAX_NUMBER_OF_OMP_THREADS
#define MAX_NUMBER_OF_OMP_THREADS 100
#endif // MAX_NUMBER_OF_OMP_THREADS

struct hrleRunTypeValues {
  static constexpr hrleSizeType SEGMENT_PT =
      std::numeric_limits<hrleSizeType>::max() -
      (MAX_NUMBER_OF_OMP_THREADS); // this value signifies end the of
                                   // the current hrleDomainSegment

  static constexpr hrleSizeType
      UNDEF_PT = // undef runs have ids up to SEGMENT_PT
      std::numeric_limits<hrleSizeType>::max() -
      (MAX_NUMBER_OF_BACKGROUND_VALUES + MAX_NUMBER_OF_OMP_THREADS);

  static constexpr hrleSizeType INACTIVE_PT =
      std::numeric_limits<hrleSizeType>::max();
};

#endif // HRLE_RUN_TYPE_VALUES_HPP
