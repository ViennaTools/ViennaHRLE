#pragma once

#include <vcVectorType.hpp>

namespace viennahrle {
using namespace viennacore;

// add different implementations for different systems if needed
typedef int IndexType;
typedef unsigned long SizeType;
typedef double CoordType;

template <int D> using Coordinate = VectorType<CoordType, D>;

template <int D> class Index {
  std::array<IndexType, D> x = {};

public:
  using value_type = IndexType;
  static constexpr int dimension = D;

  Index() = default;
  Index(const Index &v) = default;
  Index &operator=(const Index &v) = default;

  template <class T1, int D1> explicit Index(const Index<D1> &v) {
    const int k = std::min(D, D1);
    for (int i = 0; i < k; ++i)
      x[i] = v[i];
    for (int i = k; i < D; ++i)
      x[i] = 0;
  }
  explicit Index(const IndexType v[]) {
    for (int i = 0; i < D; i++)
      x[i] = v[i];
  }
  explicit Index(IndexType x0, IndexType x1, IndexType x2) {
    x[0] = x0;
    x[1] = x1;
    if constexpr (D == 3)
      x[2] = x2;
  }
  explicit Index(IndexType x0, IndexType x1) {
    x[0] = x0;
    x[1] = x1;
  }
  explicit Index(IndexType d) { fill(d); }

  template <class V> Index &operator-=(const V &v) {
    for (int i = 0; i < D; i++)
      x[i] -= v[i];
    return *this;
  }
  template <class V> Index &operator+=(const V &v) {
    for (int i = 0; i < D; i++)
      x[i] += v[i];
    return *this;
  }
  Index &operator*=(IndexType d) {
    for (int i = 0; i < D; i++)
      x[i] *= d;
    return *this;
  }
  Index &operator/=(IndexType d) {
    for (int i = 0; i < D; i++)
      x[i] /= d;
    return *this;
  }
  Index operator-() const {
    Index v;
    for (int i = 0; i < D; i++)
      v.x[i] = -x[i];
    return v;
  }

  template <class V> bool operator==(const V &v) const {
    for (int i = D - 1; i >= 0; --i)
      if (x[i] != v[i])
        return false;
    return true;
  }
  template <class V> bool operator!=(const V &v) const {
    for (int i = D - 1; i >= 0; --i)
      if (x[i] != v[i])
        return true;
    return false;
  }
  template <class V> bool operator<(const V &v) const {
    for (int i = D - 1; i >= 0; --i) {
      if (x[i] < v[i])
        return true;
      if (x[i] > v[i])
        return false;
    }
    return false;
  }
  template <class V> bool operator<=(const V &v) const { return !(*this > v); }
  template <class V> bool operator>(const V &v) const {
    for (int i = D - 1; i >= 0; --i) {
      if (x[i] > v[i])
        return true;
      if (x[i] < v[i])
        return false;
    }
    return false;
  }
  template <class V> bool operator>=(const V &v) const { return !(*this < v); }

  IndexType &operator[](int i) { return x[i]; }
  const IndexType &operator[](int i) const { return x[i]; }
  IndexType &at(int i) { return x.at(i); }
  const IndexType &at(int i) const { return x.at(i); }

  template <class V> Index &operator=(const V &v) {
    for (int i = 0; i < D; i++)
      x[i] = v[i];
    return *this;
  }

  IndexType size() const { return x.size(); }
  IndexType max_element() const {
    return *std::max_element(std::begin(x), std::end(x));
  }
  void sort() { std::sort(x.begin(), x.end()); }
  void reverse_sort() {
    std::sort(x.begin(), x.end(), std::greater<IndexType>());
  }

  IndexType *data() { return x.data(); }
  const IndexType *data() const { return x.data(); }

  const std::array<IndexType, D> &get() const { return x; }
  std::array<IndexType, D> &get() { return x; }

  auto begin() { return x.begin(); }
  auto end() { return x.end(); }
  auto begin() const { return x.begin(); }
  auto end() const { return x.end(); }

  void swap(Index &v) noexcept { x.swap(v.x); }
  void fill(IndexType value) { std::fill(x.begin(), x.end(), value); }

  std::string to_string() const {
    std::string str = "[";
    for (int i = 0; i < D - 1; i++)
      str += std::to_string(x[i]) + ", ";
    str += std::to_string(x[D - 1]) + "]";
    return str;
  }

  struct hash {
  private:
    static size_t hash_combine(size_t lhs, size_t rhs) {
      lhs ^= rhs + 0x9e3779b9 + (lhs << 6) + (lhs >> 2);
      return lhs;
    }

  public:
    size_t operator()(const Index &v) const {

      /*
        https://stackoverflow.com/questions/5889238/why-is-xor-the-default-way-to-combine-hashes
      */
      size_t result = std::hash<IndexType>{}(v[0]);
      result = hash_combine(result, std::hash<IndexType>{}(v[1]));
      if (D == 3) {
        result = hash_combine(result, std::hash<IndexType>{}(v[2]));
      }
      return result;
    }
  };
};

#define _define_operator(op)                                                   \
  /* index op index */                                                         \
  inline Index<2> operator op(const Index<2> &a, const Index<2> &b) {          \
    return Index<2>{a[0] op b[0], a[1] op b[1]};                               \
  }                                                                            \
                                                                               \
  inline Index<3> operator op(const Index<3> &a, const Index<3> &b) {          \
    return Index<3>{a[0] op b[0], a[1] op b[1], a[2] op b[2]};                 \
  }                                                                            \
                                                                               \
  /* index op scalar */                                                        \
  inline Index<2> operator op(const Index<2> &a, const IndexType & b) {        \
    return Index<2>{a[0] op b, a[1] op b};                                     \
  }                                                                            \
                                                                               \
  inline Index<3> operator op(const Index<3> &a, const IndexType & b) {        \
    return Index<3>{a[0] op b, a[1] op b, a[2] op b};                          \
  }                                                                            \
                                                                               \
  /* scalar op index */                                                        \
  inline Index<2> operator op(const IndexType & a, const Index<2> &b) {        \
    return Index<2>{a op b[0], a op b[1]};                                     \
  }                                                                            \
                                                                               \
  inline Index<3> operator op(const IndexType & a, const Index<3> &b) {        \
    return Index<3>{a op b[0], a op b[1], a op b[2]};                          \
  }

_define_operator(+);
_define_operator(-);

#undef _define_operator

template <class S, int D> S &operator<<(S &o, const Index<D> &v) {
  o << "[" << v[0];
  for (size_t i = 1; i < D; ++i)
    o << ", " << v[i];
  o << "]";
  return o;
}

template <int D> Index<D> BitMaskToIndex(unsigned int i) {
  Index<D> tmp(0);
  for (unsigned int k = 0; k < D; k++) {
    if (((1 << k) & i) != 0)
      ++tmp[k];
  }
  return tmp;
}

template <int D> int Compare(const Index<D> &v1, const Index<D> &v2) {
  for (int i = D - 1; i >= 0; --i) {
    if (v1[i] > v2[i])
      return 1;
    if (v1[i] < v2[i])
      return -1;
  }
  return 0;
}
} // namespace viennahrle
