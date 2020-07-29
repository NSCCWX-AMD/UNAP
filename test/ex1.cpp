#include "unapVector.hpp"

//- test for vectors

using namespace UNAP;

int main()
{
  Vector<scalar> v1(10);
  v1 = 1.1;
  UNAPCOUT << "OK, size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  Vector<label> v2;
  UNAPCOUT << "v2 size is " << v2.size() << ENDL;

  Vector<scalar> v3(v1);
  UNAPCOUT << "v3 size is " << v3.size() << ", 5th is " << v3[5] << ENDL;

  scalar *array1 = new scalar[12];
  forAll(i, 12) { array1[i] = (scalar)i + 0.2; }

  Vector<scalar> v4(array1, 10);
  UNAPCOUT << "v4 size is " << v4.size() << ", 5th is " << v4[5] << ENDL;

  // v1 = v4;
  // UNAPCOUT << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  // v1 += v3;
  // UNAPCOUT << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  // v1 -= 0.33;
  // UNAPCOUT << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  Vector<scalar> v5(v1 - v4);
  UNAPCOUT << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;
  UNAPCOUT << "v5 size is " << v5.size() << ", 5th is " << v5[5] << ENDL;

  scalarVector v6(10);
  v6 = v1 - v4;
  UNAPCOUT << "v6 size is " << v6.size() << ", 5th is " << v6[5] << ENDL;

  forAll(i, 10) { UNAPCOUT << "v6 " << i << ", val = " << v6[i] << ENDL; }

  scalar sum = v6.Sum();
  scalar sumMag = v6.SumMag();
  UNAPCOUT << "v6 sum = " << sum << ", sumMag = " << sumMag << ENDL;

  label newSize = 5;
  v6.SET_size(newSize);
  forAll(i, newSize)
  {
    UNAPCOUT << "after new size, v6 " << i << ", val = " << v6[i] << ENDL;
  }

  newSize = 15;
  v6.SET_size(newSize);
  forAll(i, newSize)
  {
    UNAPCOUT << "after new size, v6 " << i << ", val = " << v6[i] << ENDL;
  }

  delete[] array1;
  return 0;
}
