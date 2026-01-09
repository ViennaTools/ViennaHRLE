#ifndef HRLE_ALLOCATION_TYPE_HPP
#define HRLE_ALLOCATION_TYPE_HPP

#include <vcVectorType.hpp>

namespace viennahrle {
using namespace viennacore;
template <class SizeType = unsigned int, int D = 3> class AllocationType {
public:
  VectorType<SizeType, D> num_values;
  VectorType<SizeType, D> num_runs;

  static constexpr double allocationFactor = 1.2;

  template <class X> AllocationType &operator*=(const X &x) {
    for (int i = 0; i < D; ++i)
      num_values[i] = static_cast<SizeType>(num_values[i] * x);
    for (int i = 0; i < D; ++i)
      num_runs[i] = static_cast<SizeType>(num_runs[i] * x);
    return *this;
  }

  AllocationType &operator+=(const AllocationType &x) {
    num_values += x.num_values;
    num_runs += x.num_runs;
    return *this;
  }

  template <class X> AllocationType &operator/=(const X &x) {
    for (int i = 0; i < D; ++i)
      num_values[i] = static_cast<SizeType>(num_values[i] / x) + 1;
    for (int i = 0; i < D; ++i)
      num_runs[i] = static_cast<SizeType>(num_runs[i] / x) + 1;
    return *this;
  }

  template <class X> [[nodiscard]] AllocationType operator*(const X &x) const {
    AllocationType tmp(*this);
    tmp *= x;
    return tmp;
  }

  template <class X> [[nodiscard]] AllocationType operator/(const X &x) const {
    AllocationType tmp(*this);
    tmp /= x;
    return tmp;
  }

  template <class X> [[nodiscard]] AllocationType operator+(const X &x) const {
    AllocationType tmp(*this);
    tmp += x;
    return tmp;
  }

  AllocationType() {
    std::fill(num_runs.begin(), num_runs.end(), SizeType(0));
    std::fill(num_values.begin(), num_values.end(), SizeType(0));
  }
};
} // namespace viennahrle

#endif // HRLE_ALLOCATION_TYPE_HPP
