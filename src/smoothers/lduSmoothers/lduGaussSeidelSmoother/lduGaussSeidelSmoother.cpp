#include "lduGaussSeidelSmoother.hpp"

void UNAP::lduGaussSeidelSmoother::smooth(scalarVector &x,
                                          const lduMatrix &A,
                                          const scalarVector &b,
                                          const label nSweeps) const
{
  scalar *xPtr = x.begin();

  const label nCells = x.size();

  scalarVector bPrime(nCells, 0);

  scalar *bPrimePtr = bPrime.begin();

  const scalar *const diagPtr = A.diag().begin();

  const scalar *const upperPtr = A.upper().begin();
  const scalar *const lowerPtr = A.lower().begin();

  const label *const uPtr = A.upperAddr().begin();

  const label *const ownStartPtr = A.ownerStartAddr().begin();

  // Parallel boundary initialization.  The parallel boundary is treated
  // as an effective Jacobi interface in the boundary.
  // Note: there is a change of sign in the coupled
  // interface update.  The reason for this is that the
  // internal coefficients are all located at the l.h.s. of
  // the matrix whereas the "implicit" coefficients on the
  // coupled boundaries are all created as if the
  // coefficient contribution is of a b-kind (i.e. they
  // have a sign as if they are on the r.h.s. of the matrix.
  // To compensate for this, it is necessary to turn the
  // sign of the contribution.

  for (label sweep = 0; sweep < nSweeps; sweep++)
  {
    bPrime = 0.0;

    A.initInterfaces(x);
    A.updateInterfaces(bPrime);

    bPrime = b - bPrime;

    scalar curX;
    label fStart;
    label fEnd = ownStartPtr[0];

    forAll(cellI, nCells)
    {
      //- start and end of this row
      fStart = fEnd;
      fEnd = ownStartPtr[cellI + 1];

      //- get the accumulated neighbor side
      curX = bPrimePtr[cellI];

      //- accumulate the owner product side
      for (label curFace = fStart; curFace < fEnd; curFace++)
      {
        curX -= upperPtr[curFace] * xPtr[uPtr[curFace]];
      }

      //- finish current x
      curX /= diagPtr[cellI];

      //- distribute the neighbor side using current x
      for (label curFace = fStart; curFace < fEnd; curFace++)
      {
        bPrimePtr[uPtr[curFace]] -= lowerPtr[curFace] * curX;
      }

      xPtr[cellI] = curX;
    }
  }
}

void UNAP::lduGaussSeidelSmoother::smooth(scalarVector &x,
                                          const matrix &A,
                                          const scalarVector &b,
                                          const label nSweeps) const
{
  smooth(x, static_cast<const lduMatrix &>(A), b, nSweeps);
}
