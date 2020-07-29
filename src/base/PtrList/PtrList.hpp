/// \file PtrList.hpp
/// \brief brief information to be added
///
/// Detailed information to be added
///
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30

#ifndef PTRLIST_HPP
#define PTRLIST_HPP

#include "unap.hpp"

namespace UNAP
{
template <typename T>
class PtrList
{
private:
  label size_;

  T **ptrs_;

public:
  //- constructor
  //- null constructor
  PtrList();

  //- construct with size specified
  PtrList(const label);

  //- copy constructor
  PtrList(const PtrList<T> &);

  //- destructor
  ~PtrList()
  {
    if (ptrs_)
    {
      forAll(i, size_)
      {
        if (ptrs_[i])
        {
          delete ptrs_[i];
          ptrs_[i] = NULL;
        }
      }

      delete[] ptrs_;
      ptrs_ = NULL;
    }
  }

  //- return the number of elements in the PtrList
  inline label size() const;

  //- set size of PtrList
  inline void SET_size(const label i);

  //- return true if the PtrList is empty (ie, size() is zero)
  inline bool isEmpty() const;

  //- deep copy
  PtrList &operator=(const PtrList<T> &);

  //- return element const reference
  inline T &operator[](const label i) const { return *ptrs_[i]; }

  //- make the ptrs_[i] point to the object
  void setLevel(const label i, T &obj);

  //- delete a level
  void removeLevel(const label i);
};

template <class T>
PtrList<T>::PtrList() : size_(50), ptrs_(NULL)
{
  ptrs_ = new T *[size_];
  forAll(i, size_) { ptrs_[i] = NULL; }
}

template <class T>
PtrList<T>::PtrList(const label size) : size_(size), ptrs_(NULL)
{
  ptrs_ = new T *[size_];
  forAll(i, size_) { ptrs_[i] = NULL; }
}

template <class T>
PtrList<T>::PtrList(const PtrList<T> &oldObj) : size_(-1), ptrs_(NULL)
{
  size_ = oldObj.size();
  ptrs_ = new T *[size_];
  forAll(i, size_)
  {
    T &newObj = *ptrs_[i];
    newObj = oldObj;
  }
}

template <class T>
PtrList<T> &PtrList<T>::operator=(const PtrList<T> &oldObj)
{
  if (this == &oldObj)
  {
    return *this;
  }
  else
  {
    //- delete existed
    if (this->ptrs_)
    {
      forAll(i, this->size_)
      {
        if (this->ptrs_[i])
        {
          delete this->ptrs_[i];
          this->ptrs_[i] = NULL;
        }
      }

      delete[] this->ptrs_;
      this->ptrs_ = NULL;
    }

    //- copy
    this->size_ = oldObj.size();
    this->ptrs_ = new T *[this->size_];
    forAll(i, this->size_)
    {
      T &newObj = *(this->ptrs_[i]);
      newObj = oldObj;
    }
  }
}

template <class T>
inline label PtrList<T>::size() const
{
  return size_;
}

template <class T>
inline void PtrList<T>::SET_size(const label newSize)
{
#ifdef DEBUG
  if (newSize < 0)
  {
    UNAPCOUT << "Error in PtrList SET_size: "
             << "bad new set size " << newSize << ENDL;

    ERROR_EXIT;
  }
#endif
  label oldSize = size();

  if (newSize == 0)
  {
    this->~PtrList();
  }
  else if (oldSize == newSize)
  {
    //- nothing to do
  }
  else
  {
    T **ptrsOld = ptrs_;
    ptrs_ = new T *[newSize];

    forAll(i, newSize)
    {
      if (i < oldSize)
      {
        ptrs_[i] = ptrsOld[i];
      }
      else
      {
        ptrs_[i] = NULL;
      }
    }
    size_ = newSize;
    delete[] ptrsOld;
  }
}

template <class T>
inline bool PtrList<T>::isEmpty() const
{
  if (size_ == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

template <class T>
void PtrList<T>::removeLevel(const label leveli)
{
  label oldSize = size_;
  if (leveli < size_)
  {
    delete ptrs_[leveli];
    for (label i = leveli; i < size_ - 1; i++)
    {
      ptrs_[i] = ptrs_[i + 1];
    }
    ptrs_[size_ - 1] = NULL;
    size_--;
  }
  else
  {
    UNAPCOUT << "Error: removeLevel failed!" << ENDL;
    UNAPCOUT << "Number of levels is " << oldSize
             << ", while the level tended to delete is " << leveli << ENDL;
    ERROR_EXIT;
  }
}

template <typename T>
void PtrList<T>::setLevel(const label leveli, T &obj)
{
  if (leveli < size_)
  {
    ptrs_[leveli] = &obj;
  }
  else
  {
    UNAPCOUT << "Error: setLevel failed!" << ENDL;
    UNAPCOUT << "Number of levels is " << size_
             << ", while the level tended to set is " << leveli << ENDL;
    ERROR_EXIT;
  }
}

}  // namespace UNAP

#endif  //- PTRLIST_HPP
