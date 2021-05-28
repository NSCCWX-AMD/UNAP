#include "lduMatrix.hpp"

//- test for lduMatrix

using namespace UNAP;

int main()
{
  label nCells = 10;
  scalarVector x(nCells, 0.0);

  forAll(i, nCells)
  {
    x[i] = ((scalar)i + 0.2) / 0.4;
    std::cout << "x = " << x[i] << ENDL;
  }

  labelVector lowerAddr(nCells);
  labelVector upperAddr(nCells);

  forAll(i, nCells)
  {
    lowerAddr[i] = i == nCells - 1 ? 0 : i;
    upperAddr[i] = i == nCells - 1 ? nCells - 1 : i + 1;
  }

  scalarVector lower(nCells, -1.0);
  scalarVector upper(nCells, -1.0);
  scalarVector diag(nCells, 4.0);

  lduMatrix lduM(nCells, lowerAddr, upperAddr, lower, diag, upper);

  matrix *A = &lduM;

  scalarVector Apsi(nCells, 0.0);
  A->spMV(Apsi, x);

  forAll(i, nCells)
  {
    std::cout << "i = " << i << ", Apsi = " << Apsi[i] << ENDL;
  }

  //- for checking the correct
}