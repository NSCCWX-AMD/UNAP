#include "unapMatrix.hpp"

//- no preconditioner
void UNAP::matrix::preconditioner::precondition(scalarVector &w,
                                                const scalarVector &r) const
{
  const label nCells = w.size();
  const scalar *rPtr = r.begin();

  forAll(i, nCells) { w[i] = rPtr[i]; }
}
