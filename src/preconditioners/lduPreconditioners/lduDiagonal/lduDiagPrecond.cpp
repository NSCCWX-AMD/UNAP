#include "lduDiagPrecond.hpp"

#include "lduMatrix.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

UNAP::lduDiagPrecond::lduDiagPrecond(const lduMatrix &A) : rD_(A.nCells())
{
  const label nCells = A.nCells();
  const scalar *DPtr = A.diag().begin();

  IFNOT_SWACC
  {
    forAll(i, nCells) { rD_[i] = 1.0 / DPtr[i]; }
  }
#ifdef SW_SLAVE
  else
  {
    MVM_Arrays arrays1;
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = rD_.values();
    arrays1.A2Ptr = A.diag().values();
    // rD = 1 / D
    vectorOps_host(&arrays1, &slave_userFunc_aE1Db);
  }
#endif
}

void UNAP::lduDiagPrecond::precondition(scalarVector &w,
                                        const scalarVector &r) const
{
  const label nCells = w.size();
  const scalar *rPtr = r.begin();
  const scalar *rDPtr = rD_.begin();

  IFNOT_SWACC
  {
    forAll(i, nCells) { w[i] = rDPtr[i] * rPtr[i]; }
  }
#ifdef SW_SLAVE
  else
  {
    MVM_Arrays arrays1;
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = w.values();
    arrays1.A2Ptr = rD_.values();
    arrays1.A3Ptr = r.values();
    // wPtr = rDPtr * rPtr
    vectorOps_host(&arrays1, &slave_userFunc_aEbMuc);
  }
#endif
}
