#ifndef HRLE_VECTOR_TYPE_HPP
#define HRLE_VECTOR_TYPE_HPP

#include <algorithm>
#include <cmath>
#include <functional>

template <class T, int D> class hrleVectorType {
  T x[D];

public:
  typedef T value_type;

  hrleVectorType(const hrleVectorType<T, D> &);
  template <class T1, int D1>
  explicit hrleVectorType(const hrleVectorType<T1, D1> &);
  explicit hrleVectorType(const T[]);

  hrleVectorType();

  explicit hrleVectorType(T x0, T x1, T x2) {
    x[0] = x0;
    x[1] = x1;
    if (D == 3)
      x[2] = x2;
  }
  explicit hrleVectorType(T x0, T x1) {
    x[0] = x0;
    x[1] = x1;
  }

  explicit hrleVectorType(T);
  template <class X> explicit hrleVectorType(X);

  template <class V> hrleVectorType &operator-=(const V &);
  template <class V> hrleVectorType &operator+=(const V &);
  hrleVectorType &operator*=(T);
  hrleVectorType &operator/=(T);

  template <class V> hrleVectorType<T, D> operator-(const V &) const;
  hrleVectorType<T, D> operator-() const;

  template <class V> hrleVectorType<T, D> operator+(const V &) const;

  template <class V> bool operator==(const V &) const;
  template <class V> bool operator!=(const V &) const;
  template <class V> bool operator<(const V &) const;
  template <class V> bool operator<=(const V &) const;
  template <class V> bool operator>(const V &) const;
  template <class V> bool operator>=(const V &) const;

  hrleVectorType<T, D> operator*(T) const;

  hrleVectorType<T, D> operator/(T) const;

  T &operator[](int);
  const T &operator[](int) const;

  template <class V> hrleVectorType<T, D> &operator=(const V &);

  T size() const;
  T element_max() const;
  void sort();
  void reverse_sort();
};

template <class T, int D> inline hrleVectorType<T, D>::hrleVectorType() {}

template <class T, int D>
template <class T1, int D1>
inline hrleVectorType<T, D>::hrleVectorType(const hrleVectorType<T1, D1> &v) {
  const int k = std::min(D, D1);
  for (int i = 0; i < k; i++)
    x[i] = v[i];
  for (int i = k; i < D; i++)
    x[i] = T(0);
}

template <class T, int D>
inline hrleVectorType<T, D>::hrleVectorType(const T v[]) {
  for (int i = 0; i < D; i++)
    x[i] = v[i];
}

template <class T, int D>
template <class V>
inline hrleVectorType<T, D> &hrleVectorType<T, D>::operator=(const V &v) {
  for (int i = 0; i < D; i++)
    x[i] = v[i];
  return *this;
}

template <class T, int D>
inline hrleVectorType<T, D>::hrleVectorType(const hrleVectorType<T, D> &v) {
  for (int i = 0; i < D; i++)
    x[i] = v[i];
}

template <class T, int D> inline T &hrleVectorType<T, D>::operator[](int i) {
  return x[i];
}

template <class T, int D>
inline const T &hrleVectorType<T, D>::operator[](int i) const {
  return x[i];
}

template <class T, int D>
template <class V>
inline hrleVectorType<T, D> &hrleVectorType<T, D>::operator-=(const V &v) {
  for (int i = 0; i < D; i++)
    x[i] -= v[i];
  return *this;
}

template <class T, int D>
template <class V>
inline hrleVectorType<T, D> &hrleVectorType<T, D>::operator+=(const V &v) {
  for (int i = 0; i < D; i++)
    x[i] += v[i];
  return *this;
}

template <class T, int D>
inline hrleVectorType<T, D> &hrleVectorType<T, D>::operator*=(T d) {
  for (int i = 0; i < D; i++)
    x[i] *= d;
  return *this;
}

template <class T, int D>
inline hrleVectorType<T, D> &hrleVectorType<T, D>::operator/=(T d) {
  for (int i = 0; i < D; i++)
    x[i] /= d;
  return *this;
}

template <class T, int D>
template <class V>
inline hrleVectorType<T, D> hrleVectorType<T, D>::operator-(const V &v) const {
  hrleVectorType<T, D> w;
  for (int i = 0; i < D; i++)
    w.x[i] = x[i] - v[i];
  return w;
}

template <class T, int D>
template <class V>
inline hrleVectorType<T, D> hrleVectorType<T, D>::operator+(const V &v) const {
  hrleVectorType<T, D> w;
  for (int i = 0; i < D; i++)
    w.x[i] = x[i] + v[i];
  return w;
}

template <class T, int D>
inline hrleVectorType<T, D> hrleVectorType<T, D>::operator-() const {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v.x[i] = -x[i];
  return v;
}

template <class T, int D>
inline hrleVectorType<T, D> hrleVectorType<T, D>::operator*(T d) const {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v.x[i] = x[i] * d;
  return v;
}

template <class T, int D>
inline hrleVectorType<T, D> hrleVectorType<T, D>::operator/(T d) const {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v.x[i] = x[i] / d;
  return v;
}

template <class T, int D> inline hrleVectorType<T, D>::hrleVectorType(T d) {
  for (int i = 0; i < D; i++)
    x[i] = d;
}

template <class T, int D>
template <class X>
inline hrleVectorType<T, D>::hrleVectorType(X d) {
  for (int i = 0; i < D; i++)
    x[i] = d[i];
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator==(const V &v) const {
  for (int i = D - 1; i >= 0; i--)
    if (x[i] != v[i])
      return false;
  return true;
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator!=(const V &v) const {
  for (int i = D - 1; i >= 0; i--)
    if (x[i] != v[i])
      return true;
  return false;
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator<(const V &v) const {
  for (int i = D - 1; i >= 0; i--) {
    if (x[i] < v[i])
      return true;
    else if (x[i] > v[i])
      return false;
  }
  return false;
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator>(const V &v) const {
  for (int i = D - 1; i >= 0; i--) {
    if (x[i] > v[i])
      return true;
    else if (x[i] < v[i])
      return false;
  }
  return false;
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator<=(const V &v) const {
  return !(*this > v);
}

template <class T, int D>
template <class V>
inline bool hrleVectorType<T, D>::operator>=(const V &v) const {
  return !(*this < v);
}

template <class T, int D> inline T hrleVectorType<T, D>::size() const {
  return (sizeof(x) / sizeof(*x));
}

template <class T, int D> inline T hrleVectorType<T, D>::element_max() const {
  return *std::max_element(std::begin(x), std::end(x));
}

template <class T, int D> inline void hrleVectorType<T, D>::sort() {
  std::sort(x, x + D);
}

template <class T, int D> inline void hrleVectorType<T, D>::reverse_sort() {
  std::sort(x, x + D, std::greater<T>());
}

//###################################################################

template <class T, int D>
inline hrleVectorType<T, D> Min(const hrleVectorType<T, D> &v1,
                                const hrleVectorType<T, D> &v2) {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v[i] = std::min(v1[i], v2[i]);
  return v;
}

template <class T, int D>
inline hrleVectorType<T, D> Max(const hrleVectorType<T, D> &v1,
                                const hrleVectorType<T, D> &v2) {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v[i] = std::max(v1[i], v2[i]);
  return v;
}

template <class T, int D>
inline T Element_Max(const hrleVectorType<T, D> &v1,
                     const hrleVectorType<T, D> &v2) {
  hrleVectorType<T, D> v;
  for (int i = 0; i < D; i++)
    v[i] = std::max(v1[i], v2[i]);
  return v;
}

template <class T, int D>
inline T DotProduct(const hrleVectorType<T, D> &v1,
                    const hrleVectorType<T, D> &v2) {
  T sum(0);
  for (int i = 0; i < D; i++)
    sum += v1[i] * v2[i];
  return sum;
}

template <class T, int D>
inline hrleVectorType<T, D> operator*(T d, const hrleVectorType<T, D> &v) {
  return v * d;
}

template <class T>
inline hrleVectorType<T, 2> RotateLeft(const hrleVectorType<T, 2> &v) {
  return hrleVectorType<T, 2>(-v[1], v[0]);
}

template <class T>
inline hrleVectorType<T, 2> RotateRight(const hrleVectorType<T, 2> &v) {
  return hrleVectorType<T, 2>(v[1], -v[0]);
}

template <class T>
inline hrleVectorType<T, 3> CrossProduct(const hrleVectorType<T, 3> &v1,
                                         const hrleVectorType<T, 3> &v2) {
  return hrleVectorType<T, 3>(v1[1] * v2[2] - v1[2] * v2[1],
                              v1[2] * v2[0] - v1[0] * v2[2],
                              v1[0] * v2[1] - v1[1] * v2[0]);
}

template <class T>
inline T CrossProduct(const hrleVectorType<T, 2> &v1,
                      const hrleVectorType<T, 2> &v2) {
  return v1[0] * v2[1] - v1[1] * v2[0];
}

template <int D, class T>
inline hrleVectorType<T, D> Normalize(const hrleVectorType<T, D> &v) {
  T n = Normalize(v);
  if (n <= 0.)
    return hrleVectorType<T, D>(T(0));
  return v / n;
}

template <int D, class T> inline T Normalize(const hrleVectorType<T, D> &v) {
  T max = std::abs(v[0]);
  for (int i = 1; i < D; i++)
    max = std::max(max, std::abs(v[i]));
  if (max == 0.)
    return T(0.);
  T d = 0.;
  for (int i = 0; i < D; i++) {
    T t(v[i] / max);
    d += t * t;
  }
  return max * std::sqrt(d);
}

template <int D, class T>
inline T Norm2(
    const hrleVectorType<T, D> &v) { // squared l2 norm TODO name is misleading
  return Dot(v, v);
}

template <int D, class T>
inline T NormL2(const hrleVectorType<T, D> &v) { // l2 norm
  return std::sqrt(Norm2(v));
}

template <int D, class T>
inline int ManhattanNorm(const hrleVectorType<T, D> &v) {
  T k(0);
  for (int i = 0; i < D; i++)
    k += std::abs(v[i]);
  return k;
}

// template <class V, class T,int D> inline bool operator<(const V& v0,const
// hrleVectorType<T,D>& v1) {
//     return (v1.operator>(v0));
// }

// template <class V, class T,int D> inline bool operator>(const V& v0,const
// hrleVectorType<T,D>& v1) {
//     return (v1.operator<(v0));
// }

// template <class V, class T,int D> inline bool operator<=(const V& v0,const
// hrleVectorType<T,D>& v1) {
//     return (v1.operator>=(v0));
// }

// template <class V, class T,int D> inline bool operator>=(const V& v0,const
// hrleVectorType<T,D>& v1) {
//     return (v1.operator<=(v0));
// }

// template <class V, class T,int D> inline hrleVectorType<T,D> operator+(const
// V& v0,const hrleVectorType<T,D>& v1) {
//     return v1.operator+(v0);
// }

// template <class V, class T,int D> inline hrleVectorType<T,D> operator-(const
// V& v0,const hrleVectorType<T,D>& v1) {
//     return (-v1).operator+(v0);
// }

template <class T> inline T Volume(const hrleVectorType<T, 2> *p) {
  return ((p[1] - p[0]) % (p[2] - p[0]));
}

template <class T> inline T Volume(const hrleVectorType<T, 3> *p) {
  return ((p[1] - p[0]) * ((p[2] - p[0]) % (p[3] - p[0])));
}

/*template <class T, int D> inline T Volume(const hrleVectorType<T,D>* p) {
    return Det(mat::mat<T,D>(p));
}*/

/*template <class T> inline hrleVectorType<T,2> NormalhrleVectorTypetor(const
hrleVectorType<T,2>* v) { return RotateLeft(Normalize(v[1]-v[0]));
}

template <class T> inline hrleVectorType<T,3> NormalhrleVectorTypetor(const
hrleVectorType<T,3>* c) { //TODO

    int max_e(0);
    T max_edge(0);

    for (int e=0;e<3;e++) {
        for (int coord=0;coord<3;coord++) {
            T edge=c[(e+1)%3][coord]-c[e][coord];
            if (max_edge<edge) {
                max_edge=edge;
                max_e=e;
            }
        }
    }

    if (max_edge==T(0)) return hrleVectorType<T,3>(T(0));

    return Normalize(CrossProduct(
                    (c[max_e]      -c[(max_e+2)%3])/max_edge,
                    (c[(max_e+1)%3]-c[(max_e+2)%3])
                ));
}*/

template <int D> inline int Parity(const hrleVectorType<int, D> &v) {
  return Norm(v) % 2;
}

template <int D, class T>
inline hrleVectorType<T, D> BitMaskToVector(unsigned int i) {
  hrleVectorType<T, D> tmp(T(0));
  for (unsigned int k = 0; k < D; k++) {
    if (((1 << k) & i) != 0)
      tmp[k]++;
  }
  return tmp;
}

/*template <class T, int D> inline int hrleVectorType<T,D>::GetDirection() const
{ int dir=-1; int num0=0; for (int i=0;i<D;i++) { if (x[i]!=T(0)) { num0++;
            dir=i;
        }
    }
    if (num0==1) return dir; else return -1;
}*/

template <class T, int D> inline int MinIndex(const hrleVectorType<T, D> &v) {
  int idx = 0;
  for (int i = 1; i < D; i++) {
    if (v[i] < v[idx])
      idx = i;
  }
  return idx;
}

template <class T, int D> inline int MaxIndex(const hrleVectorType<T, D> &v) {
  int idx = 0;
  for (int i = 1; i < D; i++) {
    if (v[i] > v[idx])
      idx = i;
  }
  return idx;
}

template <class T> inline bool Orientation(const hrleVectorType<T, 3> *v) {
  return (DotProduct(CrossProduct(v[1] - v[0], v[2] - v[0]), v[3] - v[0]) >=
          -0.);
}

template <class T> inline bool Orientation(const hrleVectorType<T, 2> *v) {
  return (DotProduct(RotateLeft(v[1] - v[0]), v[2] - v[0]) >= -0.);
}

template <class S, class T, int D>
S &operator<<(S &s, const hrleVectorType<T, D> &v) {
  s << "[" << (v[0]);
  for (int i = 1; i < D; ++i)
    s << "," << (v[i]);
  s << "]";
  return s;
}

template <class T, int D>
int compare(const hrleVectorType<T, D> &v1, const hrleVectorType<T, D> &v2) {
  for (int i = D - 1; i >= 0; --i) {
    if (v1[i] > v2[i])
      return 1;
    if (v1[i] < v2[i])
      return -1;
  }
  return 0;
}

template <class T, int D>
bool AnyEqualElement(const hrleVectorType<T, D> &v1,
                     const hrleVectorType<T, D> &v2) {
  for (int i = 0; i < D - 1; ++i) {
    if (v1[i] == v2[i])
      return true;
  }
  return false;
}

/*template <class T> inline const hrleVectorType<T,3>& Normal(const
hrleVectorType<T,3>* v) { return
DotProduct(CrossProduct(v[1]-v[0],v[2]-v[0]),v[3]-v[0])
}*/

/*template <class T, int D> inline hrleVectorType<T,D>
hrleVectorType<T,D>::ReplaceNth(int i,const T& val) const { hrleVectorType<T,D>
v(x); v[i]=val; return v;
}*/

#endif // HRLE_VECTOR_TYPE_HPP
