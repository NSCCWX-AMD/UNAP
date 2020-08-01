/// \file unapVector.hpp
/// \brief brief information to be added
///
/// Detailed information to be added
///
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30

#ifndef UNAPVECTOR_HPP
#define UNAPVECTOR_HPP

#include <string.h>

#include "unap.hpp"
#include "unapMPI.hpp"

namespace UNAP
{
template <typename T>
class Vector
{
private:
  label length_;
  T *values_;

public:
  //- constructors
  //- length(0), values(NULL)
  Vector();

  Vector(const label length);

  //- copy from an existing vector
  Vector(const Vector<T> &v);

  //- copy from an existing array with a giving length
  Vector(const T *val, const label &length);

  //- copy from an existing array with a giving length
  Vector(const T *val, const label &length, const bool reUse);

  //- build a vector with a given length and same value
  Vector(const label length, const T value);

  //- destructor
  virtual ~Vector();

  //- return ith value
  inline T &operator[](const label i) { return (this->values_[i]); }

  inline const T operator[](const label i) const { return (this->values_[i]); }

  //- copy from an existing vector
  Vector &operator=(const Vector<T> &v);

  //- set all values to the scalar a
  Vector &operator=(const T &a);

  //- both vectors must be of equal length
  Vector &operator+=(const Vector<T> &v);
  Vector &operator-=(const Vector<T> &v);

  Vector &operator+=(const T &a);
  Vector &operator-=(const T &a);
  Vector &operator*=(const T &a);

  const Vector operator+(const Vector<T> &v) const;
  const Vector operator-(const Vector<T> &v) const;
  const Vector operator*(const T &a) const;

  //- return vector length
  label size() const { return (this->length_); }

  void SET_size(const label newSize);

  T *values() const { return (this->values_); }

  T Sum() const;
  T SumMag() const;
  T SumSqr() const;
  scalar SumSqrt() const;

  T *begin() const { return &(this->values_[0]); }

  T *end() const { return &(this->values_[length_ - 1]); }

  void SET_zero();
};

typedef Vector<label> labelVector;
typedef Vector<scalar> scalarVector;

scalar dot(const scalarVector &v1, const scalarVector &v2);

template <typename T>
Vector<T>::Vector() : length_(0), values_(NULL)
{
}

template <typename T>
Vector<T>::Vector(const label length) : length_(length), values_(NULL)
{
  if (length_ > 0)
  {
    values_ = new T[length_];
    forAll(i, length_) { values_[i] = 0; }
  }
}

template <typename T>
Vector<T>::Vector(const Vector<T> &v) : length_(v.length_), values_(NULL)
{
  if (length_ > 0)
  {
    values_ = new T[length_];
    memcpy(values_, v.values_, length_ * sizeof(T));
  }
}

template <typename T>
Vector<T>::Vector(const T *val, const label &length)
    : length_(length), values_(NULL)
{
  if (length_ > 0)
  {
    values_ = new T[length_];
    memcpy(values_, val, length_ * sizeof(T));
  }
}

template <typename T>
Vector<T>::Vector(const T *val, const label &length, const bool reUse)
    : length_(length), values_(NULL)
{
  if (reUse)
  {
    values_ = const_cast<T *>(val);
  }
  else
  {
    if (length_ > 0)
    {
      values_ = new T[length_];
      memcpy(values_, val, length_ * sizeof(T));
    }
  }
}

template <typename T>
Vector<T>::Vector(const label len, const T value) : length_(len), values_(NULL)
{
  if (length_ > 0)
  {
    values_ = new T[length_];
    forAll(i, length_) { values_[i] = value; }
  }
}

template <typename T>
Vector<T>::~Vector()
{
  if (values_)
  {
    delete[] values_;
    values_ = NULL;
  }
}

template <typename T>
Vector<T> &Vector<T>::operator=(const Vector<T> &v)
{
  if (&v == this)
  {
    return *this;
  }
  else
  {
    if (length_ != v.length_)
    {
      if (this->values_)
      {
        delete[] values_;
        values_ = NULL;
      }

      length_ = v.length_;
      values_ = new T[length_];
    }

    if (length_ > 0)
    {
      memcpy(values_, v.values_, length_ * sizeof(T));
    }
    return *this;
  }
}

template <typename T>
Vector<T> &Vector<T>::operator=(const T &a)
{
  if (length_ > 0)
  {
    forAll(i, length_) { values_[i] = a; }
  }
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator+=(const Vector<T> &v)
{
#ifdef DEBUG
  if (length_ != v.length_)
  {
    UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__
             << ": The length of two vectors is not same!" << ENDL;
    ERROR_EXIT;
  }
#endif

  T *val1 = values_;
  const T *val2 = v.values_;
  forAll(i, length_) { val1[i] += val2[i]; }
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator-=(const Vector<T> &v)
{
#ifdef DEBUG
  if (length_ != v.length_)
  {
    UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__
             << ": The length of two vectors is not same!" << ENDL;
    ERROR_EXIT;
  }
#endif

  T *val1 = values_;
  const T *val2 = v.values_;
  forAll(i, length_) { val1[i] -= val2[i]; }
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator+=(const T &a)
{
  if (length_ > 0)
  {
    forAll(i, length_) { values_[i] += a; }
  }
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator-=(const T &a)
{
  if (length_ > 0)
  {
    forAll(i, length_) { values_[i] -= a; }
  }
  return *this;
}

template <typename T>
Vector<T> &Vector<T>::operator*=(const T &a)
{
  if (length_ > 0)
  {
    forAll(i, length_) { values_[i] *= a; }
  }
  return *this;
}

template <typename T>
const Vector<T> Vector<T>::operator+(const Vector<T> &v) const
{
  return Vector<T>(*this) += v;
}

template <typename T>
const Vector<T> Vector<T>::operator-(const Vector<T> &v) const
{
  return Vector<T>(*this) -= v;
}

template <typename T>
const Vector<T> Vector<T>::operator*(const T &a) const
{
  return Vector<T>(*this) *= a;
}

template <typename T>
void Vector<T>::SET_size(const label newSize)
{
#ifdef DEBUG
  if (newSize < 0)
  {
    UNAPCOUT << "Error in vector SET_size: "
             << "bad new set size " << newSize << ENDL;

    ERROR_EXIT;
  }
#endif
  if (newSize == 0)
  {
    this->length_ = 0;
    if (this->values_)
    {
      delete[] this->values_;
      this->values_ = NULL;
    }
  }
  else if (newSize != this->length_)
  {
    T *nv = new T[newSize];

    if (this->length_)
    {
      label i = MIN(this->length_, newSize);
      T *vv = &this->values_[i];
      T *av = &nv[i];
      while (i--) *--av = *--vv;
    }

    if (this->values_)
    {
      delete[] this->values_;
      this->values_ = NULL;
    }

    this->length_ = newSize;
    this->values_ = nv;
  }
}

template <typename T>
T Vector<T>::Sum() const
{
  T sum = 0;
  label len = this->length_;
  T *val = this->values_;
  forAll(i, len) { sum += val[i]; }
  reduceSum(&sum);
  return sum;
}

template <typename T>
T Vector<T>::SumMag() const
{
  T sum = 0;
  label len = this->length_;
  T *val = this->values_;
  if (typeid(T) == typeid(label))
  {
    forAll(i, len) { sum += abs(val[i]); }
  }
  else if (typeid(T) == typeid(scalar))
  {
    forAll(i, len) { sum += fabs(val[i]); }
  }
  else
  {
    UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__
             << ": illegal data type!" << ENDL;
    ERROR_EXIT;
  }

  reduceSum(&sum);

  return sum;
}

template <typename T>
T Vector<T>::SumSqr() const
{
  T sum = 0;
  label len = this->length_;
  T *val = this->values_;

  forAll(i, len) { sum += val[i] * val[i]; }

  reduceSum(&sum);
  return sum;
}

template <typename T>
scalar Vector<T>::SumSqrt() const
{
  T sumSqr = SumSqr();

  return sqrt(sumSqr);
}

template <typename T>
void Vector<T>::SET_zero()
{
  label len = this->length_;
  T *val = this->values_;

  forAll(i, len) { val[i] = 0; }
}

}  // namespace UNAP
#endif  //- end UNAPVECTOR_HPP
