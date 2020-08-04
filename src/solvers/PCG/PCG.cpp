#include "PCG.hpp"

#define IFPRINT if (!MYID && ifPrint_)

UNAP::PCG::PCG(Communicator *other_comm)
    : deletePrecondPtr_(false), precondPtr_(NULL), matrix::solver(other_comm)
{
  precondPtr_ = new matrix::preconditioner(other_comm);
  deletePrecondPtr_ = true;
}

UNAP::PCG::PCG(matrix::preconditioner &precond)
    : deletePrecondPtr_(false), precondPtr_(NULL), matrix::solver(NULL)
{
  if (&precond == NULL)
  {
    commcator_->log() << "ERROR in " << __FILE__ << " " << __LINE__
                      << "The preconditioner does not exist! \n"
                      << ENDL;
    ERROR_EXIT;
  }
  else
  {
    precondPtr_ = &precond;
    this->setCommunicator(precond.getCommunicator());
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

  scalarVector pA(nCells, this->commcator_);
  scalar *pAPtr = pA.values();

  scalarVector wA(nCells, this->commcator_);
  scalar *wAPtr = wA.values();

  scalarVector rA(nCells, this->commcator_);
  scalar *rAPtr = rA.values();

  scalar wArA = GREAT;
  scalar wArAold = wArA;

  //- calculate A.psi
  A.spMV(wA, x);

  //- calculate initial residual field
  forAll(i, nCells) { rAPtr[i] = bPtr[i] - wAPtr[i]; }

#ifdef DEBUG
  //- calculate normalisation factor
  scalar normFactor = this->normFactor(b);
#endif

  solverPerf.initialResidual() = this->normFactor(rA);
  solverPerf.finalResidual() = solverPerf.initialResidual();

#ifdef DEBUG
  IFPRINT
  {
    commcator_->log() << "At nIter = ";
    std::cout.width(5);
    commcator_->log() << solverPerf.nIterations();
    commcator_->log() << ",   ini res = ";
    std::cout.width(11);
    commcator_->log() << solverPerf.initialResidual();
    commcator_->log() << ",   rel res = ";
    std::cout.width(11);
    commcator_->log() << solverPerf.initialResidual() /
                             solverPerf.initialResidual();
    commcator_->log() << ",   b norm = ";
    std::cout.width(11);
    std::cout.setf(std::ios::scientific);
    commcator_->log() << normFactor << ENDL;
  }
#endif

  do
  {
    //- store previous wArA
    wArAold = wArA;

    precondPtr_->precondition(wA, rA);

    //- update search directions
    wArA = dot(wA, rA);

    // --- Test for singularity
    if (solverPerf.checkSingularity(mag(wArA)))
    {
#ifdef DEBUG
      IFPRINT { commcator_->log() << "singularity! wArA = " << wArA << ENDL; }
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
      commcator_->log() << "At nIter = ";
      std::cout.width(5);
      commcator_->log() << solverPerf.nIterations() + 1;
      commcator_->log() << ",   fin res = ";
      std::cout.width(11);
      commcator_->log() << solverPerf.finalResidual();
      commcator_->log() << ",   rel res = ";
      std::cout.width(11);
      commcator_->log() << solverPerf.finalResidual() / normFactor << ENDL;
    }
#endif
  } while ((++solverPerf.nIterations() < maxIter_ &&
            !(solverPerf.checkConvergence(
                tolerance_, relTol_, solverPerf.nIterations(), minIter_))));

  return solverPerf;
}
