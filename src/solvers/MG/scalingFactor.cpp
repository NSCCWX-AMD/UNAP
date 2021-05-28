#include "MG.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

scalar UNAP::MGSolver::scalingFactor(scalarVector &Acf,
                                     scalarVector &field,
                                     const scalarVector &source,
                                     const matrix &A) const
{
  A.spMV(Acf, field);

  scalar sfVal[2] = {0.0, 0.0};

  scalarVector &D = A.diag();
  scalar *DPtr = D.values();
  scalar *fieldPtr = field.values();
  scalar *AcfPtr = Acf.values();
  const scalar *sourcePtr = source.values();

  label nCells = field.size();

#ifdef SWTIMER
  swTimer::startTimer("scalingFactor");
#endif
  IFNOT_SWACC_SMALL
  {
    //- while the matrix-multiply done for the scaling it is
    //- possible to perform a point-Jacobi smoothing operation cheaply
    forAll(i, nCells)
    {
      sfVal[0] += sourcePtr[i] * fieldPtr[i];
      sfVal[1] += AcfPtr[i] * fieldPtr[i];
      fieldPtr[i] += (sourcePtr[i] - AcfPtr[i]) / DPtr[i];
    }
  }
#ifdef SW_SLAVE
  else
  {
    MVM_Arrays arrays1;
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = fieldPtr;
    arrays1.A2Ptr = (scalar *)sourcePtr;
    arrays1.A3Ptr = AcfPtr;
    arrays1.A4Ptr = DPtr;
    arrays1.k1Ptr = &sfVal[0];
    arrays1.k2Ptr = &sfVal[1];
    scalingFactor_host(&arrays1);
  }
#endif
#ifdef SWTIMER
  swTimer::endTimer("scalingFactor");
#endif

  reduceSum(&sfVal[0], this->commcator_, 2);
  return sfVal[0] / (sfVal[1] + SMALL);
}
