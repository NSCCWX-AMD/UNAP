#include "cppTest1.hpp"

#include "cppFunc.hpp"

A::A(int i, double v) : a1_(i), a2_(v) {}

B::B(int i, double v) : A(i, v) {}

void B::display() const
{
  cout << "Here is B!"
       << " i = " << a1_ << ", v = " << a2_ << endl;
}

void cppFunc()
{
  B bb(2, 0.8);
  bb.display();
}

void cppfunc_() { cppFunc(); }

void show1_(long int *ptrptr, int *i, double *v, bool *symm)
{
  cout << "size of long int is " << sizeof(long int) << endl;
  cout << "size of B* is " << sizeof(B *) << endl;

  B *BPtr = new B(*i, *v);

  *(B **)ptrptr = BPtr;

  BPtr->display();

  if (*symm)
    cout << "symm is true! " << endl;
  else
    cout << "symm is false! " << endl;
}

void show2_(long int *ptrptr)
{
  B &bb = *((B *)*ptrptr);
  bb.display();
  delete (B *)*ptrptr;
}
