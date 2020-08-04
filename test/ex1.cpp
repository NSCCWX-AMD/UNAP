#include "unapVector.hpp"

//- test for vectors

using namespace UNAP;

int main()
{
  Vector<scalar> v1(10);
  v1 = 1.1;
  std::cout << "OK, size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  Vector<label> v2;
  std::cout << "v2 size is " << v2.size() << ENDL;

  Vector<scalar> v3(v1);
  std::cout << "v3 size is " << v3.size() << ", 5th is " << v3[5] << ENDL;

  scalar *array1 = new scalar[12];
  forAll(i, 12) { array1[i] = (scalar)i + 0.2; }

  Vector<scalar> v4(array1, 10);
  std::cout << "v4 size is " << v4.size() << ", 5th is " << v4[5] << ENDL;

  // v1 = v4;
  // std::cout << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  // v1 += v3;
  // std::cout << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  // v1 -= 0.33;
  // std::cout << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;

  Vector<scalar> v5(v1 - v4);
  std::cout << "v1 size is " << v1.size() << ", 5th is " << v1[5] << ENDL;
  std::cout << "v5 size is " << v5.size() << ", 5th is " << v5[5] << ENDL;

  scalarVector v6(10);
  v6 = v1 - v4;
  std::cout << "v6 size is " << v6.size() << ", 5th is " << v6[5] << ENDL;

  forAll(i, 10) { std::cout << "v6 " << i << ", val = " << v6[i] << ENDL; }

  scalar sum = v6.Sum();
  scalar sumMag = v6.SumMag();
  std::cout << "v6 sum = " << sum << ", sumMag = " << sumMag << ENDL;

  label newSize = 5;
  v6.SET_size(newSize);
  forAll(i, newSize)
  {
    std::cout << "after new size, v6 " << i << ", val = " << v6[i] << ENDL;
  }

  newSize = 15;
  v6.SET_size(newSize);
  forAll(i, newSize)
  {
    std::cout << "after new size, v6 " << i << ", val = " << v6[i] << ENDL;
  }

  delete[] array1;
  return 0;
}
