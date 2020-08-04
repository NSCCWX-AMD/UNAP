#include "eigenDiagPCG.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

UNAP::eigenDiagPCG::eigenDiagPCG(const matrix &A,
                                 scalarVector &x,
                                 const scalarVector &b,
                                 const matrix::preconditioner &precond,
                                 const label nDiagPCGs)
    : maxEigenValue_(0.0)
{
  if (A.getCommunicator() != b.getCommunicator() &&
      b.getCommunicator() != x.getCommunicator())
  {
    A.getCommunicator()->log()
        << "Error" << __FILE__ << " " << __LINE__
        << "The communicators between A, b and c are different\n";
    ERROR_EXIT;
  }
  Communicator *comm = x.getCommunicator();
  scalarVector alphas(nDiagPCGs, 0.0, comm);
  scalarVector betas(nDiagPCGs, 0.0, comm);

  //- do some PCG loops to recode alphas and betas
  //- x is renewed during this process
  diagPCGLoops(A, x, b, precond, nDiagPCGs, alphas, betas);

  computeMaxEig(alphas, betas, nDiagPCGs, 1);
}

void UNAP::eigenDiagPCG::diagPCGLoops(const matrix &A,
                                      scalarVector &x,
                                      const scalarVector &b,
                                      const matrix::preconditioner &precond,
                                      const label nDiagPCGs,
                                      scalarVector &alphas,
                                      scalarVector &betas) const
{
  if (A.getCommunicator() != b.getCommunicator() &&
      b.getCommunicator() != x.getCommunicator())
  {
    A.getCommunicator()->log()
        << "Error" << __FILE__ << " " << __LINE__
        << "The communicators between A, b and c are different\n";
    ERROR_EXIT;
  }

  Communicator *comm = x.getCommunicator();

  label nCells = x.size();
  scalar *xPtr = x.begin();

  scalarVector pA(nCells, comm);
  scalar *pAPtr = pA.begin();

  scalarVector wA(nCells, comm);
  scalar *wAPtr = wA.begin();

  const scalar *rDPtr = precond.rD().begin();

  scalarVector rA(nCells, comm);
  scalar *rAPtr = rA.begin();

  const scalar *bPtr = b.begin();

  scalar wArA = GREAT;
  scalar wArAold = wArA;

  //- calculate A.psi
  A.spMV(wA, x);

#ifdef SW_SLAVE
  MVM_Arrays arrays1;
#endif

  //- calculate initial residual field
  IFNOT_SWACC
  {
    forAll(cell, nCells) { rAPtr[cell] = bPtr[cell] - wAPtr[cell]; }
  }
#ifdef SW_SLAVE
  else
  {
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = rAPtr;
    arrays1.A2Ptr = (scalar *)bPtr;
    arrays1.A3Ptr = wAPtr;
    vectorOps_host(&arrays1, &slave_userFunc_aEbMic);
  }
#endif

  forAll(nIter, nDiagPCGs)
  {
    //- store previous wArA
    wArAold = wArA;
    wArA = 0.0;
    IFNOT_SWACC
    {
      forAll(cell, nCells)
      {
        wAPtr[cell] = rDPtr[cell] * rAPtr[cell];
        wArA += wAPtr[cell] * rAPtr[cell];
      }
    }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A1Ptr = wAPtr;
      arrays1.A2Ptr = (scalar *)rDPtr;
      arrays1.A3Ptr = rAPtr;
      arrays1.k1Ptr = &wArA;
      gSum_host(&arrays1, &slave_userFunc_digPrecondSum);
    };
#endif
    reduceSum(&wArA, comm);

    if ((mag(wArA)) < VSMALL)
    {
#ifdef DEBUG

      comm->log() << "In " << __FILE__ << " " << __LINE__ << ENDL;
      comm->log() << "Warning: singularity in calculating eigenvalue! " << ENDL;
      comm->log() << "The value of search directions is too small:  wArA = "
                  << wArA << ENDL;

#endif
      break;
    }

    if (nIter == 0)
    {
      betas[nIter] = 0.0;

      IFNOT_SWACC
      {
        forAll(cell, nCells) { pAPtr[cell] = wAPtr[cell]; }
      }
#ifdef SW_SLAVE
      else
      {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = pAPtr;
        arrays1.A2Ptr = wAPtr;
        vectorCopy_host(&arrays1);
      }
#endif
    }
    else
    {
      scalar beta = wArA / wArAold;
      betas[nIter] = beta;

      IFNOT_SWACC
      {
        forAll(cell, nCells) { pAPtr[cell] = wAPtr[cell] + beta * pAPtr[cell]; }
      }
#ifdef SW_SLAVE
      else
      {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = pAPtr;
        arrays1.A2Ptr = wAPtr;
        arrays1.k1 = beta;
        vectorOps_host(&arrays1, &slave_userFunc_aEbPk1Mua);
      }
#endif
    }

    //- update preconditioned residual
    A.spMV(wA, pA);

    scalar wApA = 0.0;

    IFNOT_SWACC
    {
      forAll(cell, nCells) { wApA += wAPtr[cell] * pAPtr[cell]; }
    }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A2Ptr = wAPtr;
      arrays1.A3Ptr = pAPtr;
      arrays1.k1Ptr = &wApA;
      gSum_host(&arrays1, &slave_userFunc_sumProd);
    }
#endif
    reduceSum(&wApA, comm);

    //- update solution and residual
    scalar alpha = wArA / wApA;
    alphas[nIter] = alpha;

    IFNOT_SWACC
    {
      forAll(cell, nCells)
      {
        xPtr[cell] += alpha * pAPtr[cell];
        rAPtr[cell] -= alpha * wAPtr[cell];
      }
    }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A1Ptr = xPtr;
      arrays1.A2Ptr = rAPtr;
      arrays1.A3Ptr = pAPtr;
      arrays1.A4Ptr = wAPtr;
      arrays1.k1 = alpha;
      arrays1.returnA2 = 1;
      vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Muc_bEbMik1Mud);
    }
#endif
  }
}

void UNAP::eigenDiagPCG::computeMaxEig(const scalarVector &alphas,
                                       const scalarVector &betas,
                                       const label nPCGs,
                                       const label k) const
{
  scalar xBegin, xEnd, xn;
  label s;
  const scalar SMALLBAND = 0.00001;
  std::vector<scalar> pvector(nPCGs);

  scalar **TriMatrix = allocateSym2D<scalar>(nPCGs);
  computeValueForMatrix(alphas, betas, TriMatrix, nPCGs);
  determineEigRange(TriMatrix, xBegin, xEnd, nPCGs);

  do
  {
    xn = (xBegin + xEnd) * 0.5;
    h14Sturm(TriMatrix, xn, pvector, s, nPCGs);
    pvector.insert(pvector.begin(), 1, 1);  // because p[0] = 1

    if (s >= k)
      xBegin = xn;
    else
      xEnd = xn;
  } while ((xEnd - xBegin) >= SMALLBAND);

  deleteSym2D(TriMatrix, nPCGs);

  maxEigenValue_ = xn;
}

void UNAP::eigenDiagPCG::computeValueForMatrix(const scalarVector &alphas,
                                               const scalarVector &betas,
                                               scalar **TriMatrix,
                                               const label nPCGs) const
{
  for (label i = 0; i < nPCGs; i++)
    for (label j = 0; j < nPCGs; j++)
    {
      TriMatrix[i][j] = 0;
    }

  for (label i = 0; i < nPCGs; i++)
  {
    // diagonal values
    if (i == 0)
      TriMatrix[i][i] = 1.0 / alphas[i];
    else
      TriMatrix[i][i] = 1.0 / alphas[i] + betas[i] / alphas[i - 1];

    // off-diag values
    if (i < nPCGs - 1) TriMatrix[i][i + 1] = sqrt(betas[i + 1]) / alphas[i];
    if (i > 0) TriMatrix[i][i - 1] = TriMatrix[i - 1][i];
  }
}

void UNAP::eigenDiagPCG::determineEigRange(scalar **TriMatrix,
                                           scalar &xBegin,
                                           scalar &xEnd,
                                           const label nPCGs) const
{
  scalar Tright, Tleft;
  scalar lambMax, lambMin;

  lambMin = lambMax = xBegin = xEnd = 0.0;
  for (label i = 0; i < nPCGs; i++)
  {
    Tright = (i == (nPCGs - 1)) ? 0 : TriMatrix[i][i + 1];
    Tleft = (i == 0) ? 0 : TriMatrix[i][i - 1];
    lambMax = TriMatrix[i][i] + fabs(Tright) + fabs(Tleft);
    lambMin = TriMatrix[i][i] - fabs(Tright) - fabs(Tleft);
    xBegin = (xBegin >= lambMin) ? lambMin : xBegin;
    xEnd = (xEnd <= lambMax) ? lambMax : xEnd;
  }
}

void UNAP::eigenDiagPCG::h14Sturm(scalar **TriMatrix,
                                  const scalar lamb,
                                  std::vector<scalar> &p,
                                  label &s,
                                  const label nPCGs) const
{
  label k;
  s = 0;

  p[0] = lamb - TriMatrix[0][0];
  if (p[0] < 0) s = 1;

  p[1] = (lamb - TriMatrix[1][1]) * p[0] - pow(TriMatrix[0][1], 2);
  if (p[1] * p[0] < 0) s++;

  for (k = 2; k < nPCGs; k++)
  {
    p[k] = (lamb - TriMatrix[k][k]) * p[k - 1] -
           pow(TriMatrix[k - 1][k], 2) * p[k - 2];
    if (p[k] * p[k - 1] < 0) s++;
  }
}
