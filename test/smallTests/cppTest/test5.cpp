#include <stdio.h>
#include <stdlib.h>

#include <iostream>

#include "string.h"

//- test for extern "C"

using namespace std;

class A
{
private:
public:
  void amg() { cout << "ok 1" << endl; }
  void amg(int i) { cout << "ok 2 and " << i << endl; }
  void amg(double v) { cout << "ok 3 and " << v << endl; }
};

void bgk(int i) { cout << "ok 1.1 and " << i << endl; }

void bgk(double v) { cout << "ok 3.3 and " << v << endl; }

#ifdef __cplusplus
extern "C"
{
#endif

  class B
  {
  private:
  public:
    void amg() { cout << "ok 1" << endl; }
    void amg(int i) { cout << "ok 2 and " << i << endl; }
    void amg(double v) { cout << "ok 3 and " << v << endl; }
  };

  void cgk(int i) { cout << "ok 1.1 and " << i << endl; }

  void cgk(double v) { cout << "ok 3.3 and " << v << endl; }

#ifdef __cplusplus
}
#endif

void main()
{
  int i = 2;
  double v = 0.2;
  A aa;
  aa.amg(i);
  B bb;
  bb.amg(v);
}
