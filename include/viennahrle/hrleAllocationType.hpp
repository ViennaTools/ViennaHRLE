#ifndef HRLE_ALLOCATION_TYPE_HPP
#define HRLE_ALLOCATION_TYPE_HPP

#include "hrleVectorType.hpp"

template <class SizeType = unsigned int, int D = 3> class hrleAllocationType {
public:
  hrleVectorType<SizeType, D> num_values;
  hrleVectorType<SizeType, D> num_runs;

  static constexpr double allocationFactor = 1.2;

  template <class X> hrleAllocationType &operator*=(const X &x) {
    for (int i = 0; i < D; ++i)
      num_values[i] = static_cast<SizeType>(num_values[i] * x);
    for (int i = 0; i < D; ++i)
      num_runs[i] = static_cast<SizeType>(num_runs[i] * x);
    return *this;
  }

  hrleAllocationType &operator+=(const hrleAllocationType &x) {
    num_values += x.num_values;
    num_runs += x.num_runs;
    return *this;
  }

  template <class X> hrleAllocationType &operator/=(const X &x) {
    for (int i = 0; i < D; ++i)
      num_values[i] = static_cast<SizeType>(num_values[i] / x) + 1;
    for (int i = 0; i < D; ++i)
      num_runs[i] = static_cast<SizeType>(num_runs[i] / x) + 1;
    return *this;
  }

  template <class X> hrleAllocationType operator*(const X &x) {
    hrleAllocationType tmp(*this);
    tmp *= x;
    return tmp;
  }

  template <class X> hrleAllocationType operator/(const X &x) {
    hrleAllocationType tmp(*this);
    tmp /= x;
    return tmp;
  }

  template <class X> hrleAllocationType operator+(const X &x) {
    hrleAllocationType tmp(*this);
    tmp += x;
    return tmp;
  }

  hrleAllocationType() : num_values(SizeType(0)), num_runs(SizeType(0)) {}
};

#endif // HRLE_ALLOCATION_TYPE_HPP
