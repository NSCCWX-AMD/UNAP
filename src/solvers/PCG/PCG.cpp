#include "PCG.hpp"

#define IFPRINT if (!MYID && ifPrint_)

UNAP::PCG::PCG() : deletePrecondPtr_(false), precondPtr_(NULL)
{
  precondPtr_ = new matrix::preconditioner;
  deletePrecondPtr_ = true;
}

UNAP::PCG::PCG(matrix::preconditioner &precond)
    : deletePrecondPtr_(false), precondPtr_(NULL)
{
  if (&precond == NULL)
  {
    precondPtr_ = new matrix::preconditioner;
    deletePrecondPtr_ = true;
  }
  else
  {
    precondPtr_ = &precond;
  }
}

UNAP::matrix::solverPerformance UNAP::PCG::solve(scalarVector &x,
                                                 const matrix &A,
                                                 const scalarVector &b) const
{
  matrix::solverPerformance solverPerf;

  label nCells = x.size();
  scalar *xPtr = x.values();
  const scalar *bPtr = b.values();

  scalarVector pA(nCells);
  scalar *pAPtr = pA.values();

  scalarVector wA(nCells);
  scalar *wAPtr = wA.values();

  scalarVector rA(nCells);
  scalar *rAPtr = rA.values();

  scalar wArA = GREAT;
  scalar wArAold = wArA;

  //- calculate A.psi
  A.spMV(wA, x);

  //- calculate initial residual field
  forAll(i, nCells) { rAPtr[i] = bPtr[i] - wAPtr[i]; }


  //- calculate normalisation factor
  scalar normFactor = this->normFactor(b);

  solverPerf.initialResidual() = this->normFactor(rA);
  solverPerf.finalResidual() = solverPerf.initialResidual();

#ifdef DEBUG
  IFPRINT
  {
    UNAPCOUT << "At nIter = ";
    std::cout.width(5);
    UNAPCOUT << solverPerf.nIterations();
    UNAPCOUT << ",   ini res = ";
    std::cout.width(11);
    UNAPCOUT << solverPerf.initialResidual();
    UNAPCOUT << ",   rel res = ";
    std::cout.width(11);
    UNAPCOUT << solverPerf.initialResidual() / normFactor;
    UNAPCOUT << ",   b norm = ";
    std::cout.width(11);
    std::cout.setf(std::ios::scientific);
    UNAPCOUT << normFactor << ENDL;
  }
#endif

  solverPerf.initialResidual() = normFactor;

  do
  {
    //- store previous wArA
    wArAold = wArA;

    precondPtr_->precondition(wA, rA);

    //- update search directions
    wArA = dot(wA, rA);

    // --- Test for singularity
    if (solverPerf.checkSingularity(mag(wArA) / normFactor))
    {
#ifdef DEBUG
      IFPRINT { UNAPCOUT << "singularity! wArA = " << wArA << ENDL; }
#endif
      break;
    }

    if (solverPerf.nIterations() == 0)
    {
      forAll(i, nCells) { pAPtr[i] = wAPtr[i]; }
    }
    else
    {
      scalar beta = wArA / wArAold;

      forAll(i, nCells) { pAPtr[i] = wAPtr[i] + beta * pAPtr[i]; }
    }

    //- update preconditioned residual
    A.spMV(wA, pA);

    scalar wApA = dot(wA, pA);

    //- update solution and residual
    scalar alpha = wArA / wApA;

    forAll(i, nCells)
    {
      xPtr[i] += alpha * pAPtr[i];
      rAPtr[i] -= alpha * wAPtr[i];
    }

    // solverPerf.finalResidual() = rA.SumMag() / normFactor;
    solverPerf.finalResidual() = this->normFactor(rA);

#ifdef DEBUG
    IFPRINT
    {
      UNAPCOUT << "At nIter = ";
      std::cout.width(5);
      UNAPCOUT << solverPerf.nIterations() + 1;
      UNAPCOUT << ",   fin res = ";
      std::cout.width(11);
      UNAPCOUT << solverPerf.finalResidual();
      UNAPCOUT << ",   rel res = ";
      std::cout.width(11);
      UNAPCOUT << solverPerf.finalResidual() / normFactor << ENDL;
    }
#endif
  } while ((++solverPerf.nIterations() < maxIter_ &&
            !(solverPerf.checkConvergence(
                tolerance_, relTol_, solverPerf.nIterations(), minIter_))));

  return solverPerf;
}
