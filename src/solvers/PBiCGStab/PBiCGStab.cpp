#include "PBiCGStab.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

#include "printUNAP.hpp"

#define IFPRINT if (!MYID && ifPrint_)

UNAP::PBiCGStab::PBiCGStab(Communicator *other_comm)
    : deletePrecondPtr_(false), precondPtr_(NULL), matrix::solver(other_comm)
{
  precondPtr_ = new matrix::preconditioner(other_comm);
  deletePrecondPtr_ = true;
}

UNAP::PBiCGStab::PBiCGStab(matrix::preconditioner &precond)
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

UNAP::matrix::solverPerformance UNAP::PBiCGStab::solve(
    scalarVector &x, const matrix &A, const scalarVector &b) const
{
  matrix::solverPerformance solverPerf;

  if (x.getCommunicator() != A.getCommunicator() &&
      A.getCommunicator() != this->commcator_)
  {
    commcator_->log()
        << "ERROR in " << __FILE__ << " " << __LINE__
        << "The communicators between A, x and  solver are different\n"
        << ENDL;
    ERROR_EXIT;
  }
  label nCells = x.size();

  scalar *xPtr = x.values();
  const scalar *bPtr = b.values();

  scalarVector pA(nCells, this->commcator_);
  scalar *pAPtr = pA.values();

  scalarVector yA(nCells, this->commcator_);
  scalar *yAPtr = yA.values();

  scalarVector rA(nCells, this->commcator_);
  scalar *rAPtr = rA.values();

  scalarVector rA0(nCells, this->commcator_);
  scalar *rA0Ptr = rA0.values();

  //- calculate A.psi
  A.spMV(yA, x);

#ifdef SW_SLAVE
  MVM_Arrays arrays1;
#endif

  //- calculate initial residual field and store
  //- calculate norm-factor
#ifdef SWTIMER
  swTimer::startTimer("pbicgstab");
#endif
  IFNOT_SWACC
  {
    forAll(i, nCells)
    {
      rAPtr[i] = bPtr[i] - yAPtr[i];
      rA0Ptr[i] = rAPtr[i];
    }
    solverPerf.initialResidual() = this->normFactor(rA);
  }
#ifdef SW_SLAVE
  else
  {
    scalar rASqr = 0.0;
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = rAPtr;
    arrays1.A2Ptr = (scalar *)bPtr;
    arrays1.A3Ptr = yAPtr;
    arrays1.A4Ptr = rA0Ptr;
    arrays1.k1Ptr = &rASqr;
    residualNormFactor_host(&arrays1);
    reduceSum(&rASqr);
    solverPerf.initialResidual() = sqrt(rASqr);
  }
#endif
#if (SWTIMER)
  swTimer::endTimer("pbicgstab");
#endif

  // {
  //     x = rA;
  //     return solverPerf;
  // }

  solverPerf.finalResidual() = solverPerf.initialResidual();
  solverPerf.previousResidual() = solverPerf.finalResidual();

#ifdef DEBUG
  //- calculate normalisation factor
  scalar normFactor = this->normFactor(b);
  IFPRINT
  {
    commcator_->log() << "At nIter = ";
    std::cout.width(5);
    commcator_->log() << solverPerf.nIterations();
    commcator_->log() << ",   ini res = ";
    std::cout.width(11);
    std::cout.setf(std::ios::scientific);
    commcator_->log() << solverPerf.initialResidual();
    commcator_->log() << ",   rel res = ";
    commcator_->log() << solverPerf.initialResidual() /
                             solverPerf.initialResidual();
    commcator_->log() << ",   rhs  norm = ";
    commcator_->log() << normFactor << ENDL;
  }
#endif

  scalarVector AyA(nCells, this->commcator_);
  scalar *AyAPtr = AyA.values();

  scalarVector sA(nCells, this->commcator_);
  scalar *sAPtr = sA.values();

  scalarVector zA(nCells, this->commcator_);
  scalar *zAPtr = zA.values();

  scalarVector tA(nCells, this->commcator_);
  scalar *tAPtr = tA.values();

  //- initial values not used
  scalar rA0rA = 0;
  scalar alpha = 0;
  scalar omega = 0;

  do
  {
    //- store previous rA0rA
    const scalar rA0rAold = rA0rA;

    //- update search directions
#ifdef SWTIMER
    swTimer::startTimer("pbicgstab");
#endif
    IFNOT_SWACC { rA0rA = dot(rA0, rA); }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A2Ptr = rA0Ptr;
      arrays1.A3Ptr = rAPtr;
      arrays1.k1Ptr = &rA0rA;
      // rA0rA += rA0 * rA
      gSum_host(&arrays1, &slave_userFunc_sumProd);
      reduceSum(&rA0rA);
    }
#endif
#ifdef SWTIMER
    swTimer::endTimer("pbicgstab");
#endif

    // --- Test for singularity
    if (solverPerf.checkSingularity(mag(rA0rA)))
    {
#ifdef DEBUG
      IFPRINT { commcator_->log() << "singularity! rA0rA = " << rA0rA << ENDL; }
#endif
      break;
    }

    // --- update pA
    if (solverPerf.nIterations() == 0)
    {
      IFNOT_SWACC
      {
        forAll(i, nCells) { pAPtr[i] = rAPtr[i]; }
      }
#ifdef SW_SLAVE
      else
      {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = pAPtr;
        arrays1.A2Ptr = rAPtr;
        // pA = rA
        vectorCopy_host(&arrays1);
      }
#endif
    }
    else
    {
      const scalar beta = (rA0rA / rA0rAold) * (alpha / omega);
      IFNOT_SWACC
      {
        forAll(i, nCells)
        {
          pAPtr[i] = rAPtr[i] + beta * (pAPtr[i] - omega * AyAPtr[i]);
        }
      }
#ifdef SW_SLAVE
      else
      {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = pAPtr;
        arrays1.A2Ptr = rAPtr;
        arrays1.A3Ptr = AyAPtr;
        arrays1.k1 = beta;
        arrays1.k2 = omega;
        // pA = rA + beta*(pA - omega*AyA)
        vectorOps_host(&arrays1, &slave_userFunc_aEbPk1MuSaMik2MucS);
      }
#endif
    }

    //- reserved for preconditioners
    precondPtr_->precondition(yA, pA);

    //- calculate AyA
    A.spMV(AyA, yA);

    scalar rA0AyA = 0.0;

    IFNOT_SWACC { rA0AyA = dot(rA0, AyA); }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A2Ptr = rA0Ptr;
      arrays1.A3Ptr = AyAPtr;
      arrays1.k1Ptr = &rA0AyA;
      // rA0AyA += rA0 * AyA
      gSum_host(&arrays1, &slave_userFunc_sumProd);
      reduceSum(&rA0AyA);
    }
#endif

    alpha = rA0rA / (rA0AyA);

    //- calculate sA
    //- test sA for convergence

    IFNOT_SWACC
    {
      forAll(i, nCells) { sAPtr[i] = rAPtr[i] - alpha * AyAPtr[i]; }
      solverPerf.finalResidual() = this->normFactor(sA);
    }
#ifdef SW_SLAVE
    else
    {
      scalar sATemp;
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A1Ptr = sAPtr;
      arrays1.A2Ptr = rAPtr;
      arrays1.A3Ptr = AyAPtr;
      arrays1.k1Ptr = &sATemp;
      arrays1.k1 = alpha;
      gSum_host(&arrays1, &slave_userFunc_residualSumKSqr);
      reduceSum(&sATemp);
      solverPerf.finalResidual() = sqrt(sATemp);
    }
#endif

    if (solverPerf.checkConvergence(
            tolerance_, relTol_, solverPerf.nIterations() + 1, minIter_))
    {
      IFNOT_SWACC
      {
        forAll(i, nCells) { xPtr[i] += alpha * yAPtr[i]; }
      }
#ifdef SW_SLAVE
      else
      {
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = xPtr;
        arrays1.A2Ptr = yAPtr;
        arrays1.k1 = alpha;
        // x += k1 * yA
        vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Mub);
      }
#endif

      solverPerf.nIterations()++;

#ifdef DEBUG
      scalar convergenceRate =
          solverPerf.finalResidual() / solverPerf.previousResidual();
      solverPerf.previousResidual() = solverPerf.finalResidual();
      IFPRINT
      {
        commcator_->log() << "At nIter = ";
        std::cout.width(5);
        commcator_->log() << solverPerf.nIterations();
        commcator_->log() << ",   fin res = ";
        std::cout.width(11);
        std::cout.setf(std::ios::scientific);
        commcator_->log() << solverPerf.finalResidual();
        commcator_->log() << ",   rel res = ";
        commcator_->log() << solverPerf.finalResidual() / normFactor;
        commcator_->log() << ",   conv rate = ";
        commcator_->log() << convergenceRate << ENDL;
      }
#endif
      return solverPerf;
    }

    //- reserved for preconditioners
    precondPtr_->precondition(zA, sA);

    //- calculate tA
    A.spMV(tA, zA);

    scalar tAtA = 0.0;

    //- calculate omega from tA and sA
    //- (cheaper than using zA with preconditioned tA)
    IFNOT_SWACC
    {
      tAtA = tA.SumSqr();
      omega = dot(tA, sA);
    }
#ifdef SW_SLAVE
    else
    {
      scalar temp[2];
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A2Ptr = tAPtr;
      arrays1.A3Ptr = sAPtr;
      // arrays1.k1Ptr = &tAtA;
      // arrays1.k2Ptr = &omega;
      arrays1.k1Ptr = &temp[0];
      arrays1.k2Ptr = &temp[1];
      gSum_host(&arrays1, &slave_userFunc_sumSqrDot);
      // reduceSum(&temp[0], 2);
      tAtA = temp[0];
      omega = temp[1];
      reduce(tAtA);
      reduce(omega);
    }
#endif

    omega /= tAtA;
    //- update solution and residual
    IFNOT_SWACC
    {
      forAll(i, nCells)
      {
        xPtr[i] += alpha * yAPtr[i] + omega * zAPtr[i];
        rAPtr[i] = sAPtr[i] - omega * tAPtr[i];
      }

      solverPerf.finalResidual() = this->normFactor(rA);
    }
#ifdef SW_SLAVE
    else
    {
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A1Ptr = xPtr;
      arrays1.A2Ptr = yAPtr;
      arrays1.A3Ptr = zAPtr;
      arrays1.k1 = alpha;
      arrays1.k2 = omega;
      // x += k1 * yA + k2 * zA
      vectorOps_host(&arrays1, &slave_userFunc_aEaPk1MubPk2Muc);

      scalar rATemp;
      init_MVM_Arrays(&arrays1, nCells);
      arrays1.A1Ptr = rAPtr;
      arrays1.A2Ptr = sAPtr;
      arrays1.A3Ptr = tAPtr;
      arrays1.k1Ptr = &rATemp;
      arrays1.k1 = omega;
      // rA = sA - k1 * tA
      gSum_host(&arrays1, &slave_userFunc_residualSumKSqr);
      reduceSum(&rATemp);
      solverPerf.finalResidual() = sqrt(rATemp);
    }
#endif

#ifdef DEBUG
    scalar convergenceRate =
        solverPerf.finalResidual() / solverPerf.previousResidual();
    solverPerf.previousResidual() = solverPerf.finalResidual();
    IFPRINT
    {
      commcator_->log() << "At nIter = ";
      std::cout.width(5);
      commcator_->log() << solverPerf.nIterations() + 1;
      commcator_->log() << ",   fin res = ";
      std::cout.width(11);
      std::cout.setf(std::ios::scientific);
      commcator_->log() << solverPerf.finalResidual();
      commcator_->log() << ",   rel res = ";
      commcator_->log() << solverPerf.finalResidual() / normFactor;
      commcator_->log() << ",   conv rate = ";
      commcator_->log() << convergenceRate << ENDL;
    }
#endif

  } while ((++solverPerf.nIterations() < maxIter_ &&
            !(solverPerf.checkConvergence(
                tolerance_, relTol_, solverPerf.nIterations(), minIter_))));

  return solverPerf;
}
