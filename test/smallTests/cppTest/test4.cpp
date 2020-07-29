#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "string.h"

//- test for memory

using namespace std;

class A
{
private:
  int size_;
  double *ptrs_;

public:
  A(int size);
  ~A()
  {
    cout << "here is destructor from A 0." << endl;
    if (ptrs_)
    {
      delete ptrs_;
      ptrs_ = NULL;
    }
    cout << "here is destructor from A 1." << endl;
  }

  int size() const { return size_; }

  double &values(int i) { return ptrs_[i]; }
};

A::A(int size) : size_(size), ptrs_(NULL)
{
  if (size_ > 0)
  {
    ptrs_ = new double[size_];
  }
}

class B
{
private:
  A *APtr_;
  int *bPtr_;

public:
  B(A &a) : APtr_(NULL), bPtr_(NULL)
  {
    int size = a.size();
    APtr_ = new A(size);
    A &newA = *APtr_;
    for (int i = 0; i < size; ++i)
    {
      newA.values(i) = a.values(i);
    }
  }

  B(int size) : APtr_(NULL), bPtr_(NULL)
  {
    if (size > 0)
    {
      bPtr_ = new int[size];
    }
  }

  ~B()
  {
    cout << "here is destructor from B 0." << endl;
    if (APtr_)
    {
      delete APtr_;
      APtr_ = NULL;
    }

    if (bPtr_)
    {
      delete bPtr_;
      bPtr_ = NULL;
    }
    cout << "here is destructor from B 1." << endl;
  }
};

void main()
{
  A aa(10);
  B bb1(aa);
  B bb2(10);
}
