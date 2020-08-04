#include "MG.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

#include <sstream>

#include "printUNAP.hpp"

#define IFPRINT if (!MYID && ifPrint_)

UNAP::matrix::solverPerformance UNAP::MGSolver::solve(
    scalarVector &x, const matrix &A, const scalarVector &b) const
{
  //- setup class containing solver performance data
  matrix::solverPerformance solverPerf;

  const label nCells = x.size();

  //- calculate A.psi used to calculate the initial residual
  scalarVector Apsi(nCells, this->commcator_);
  A.spMV(Apsi, x);

  //- create the storage for the finestCorrection which may be used as a
  //  temporary in normFactor
  scalarVector finestCorrection(nCells, this->commcator_);

#ifdef SW_SLAVE
  MVM_Arrays arrays1;
#endif

  //- calculate initial finest-grid residual field
  //- calculate normalised residual for convergence test
  scalarVector finestResidual(nCells, this->commcator_);
  IFNOT_SWACC
  {
    finestResidual = b - Apsi;
    solverPerf.initialResidual() = this->normFactor(finestResidual);
  }
#ifdef SW_SLAVE
  else
  {
    scalar temp = 0.0;
    init_MVM_Arrays(&arrays1, nCells);
    arrays1.A1Ptr = finestResidual.values();
    arrays1.A2Ptr = b.values();
    arrays1.A3Ptr = Apsi.values();
    arrays1.k1Ptr = &temp;
    residualNormFactor_host(&arrays1);
    reduceSum(&temp);
    solverPerf.initialResidual() = sqrt(temp);
  }
#endif

  solverPerf.finalResidual() = solverPerf.initialResidual();
  solverPerf.previousResidual() = solverPerf.finalResidual();

#ifdef DEBUG
  //- calculate normalization factor
  scalar normFactor = this->normFactor(b);
  IFPRINT
  {
    commcator_->log() << "At cycle = ";
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
// swTimer::startTimer("MG Vcycle");
#endif

  if (!solverPerf.checkConvergence(
          tolerance_, relTol_, solverPerf.nIterations(), minIter_))
  {
    label coarseLevels = agglomeration_.size();

    //- create coarse grid correction fields
    PtrList<scalarVector> coarseCorrFields(coarseLevels);

    //- create coarse grid sources
    PtrList<scalarVector> coarseSources(coarseLevels);

    initVcycle(coarseCorrFields, coarseSources);

    forAll(levelI, coarseLevels) { agglomeration_.agglomerateMatrix(levelI); }

    do
    {
      Vcycle(x,
             b,
             Apsi,
             finestCorrection,
             finestResidual,
             coarseCorrFields,
             coarseSources);

      //- calculate finest level residual field
      A.spMV(Apsi, x);

      IFNOT_SWACC
      {
        finestResidual = b - Apsi;
        solverPerf.finalResidual() = this->normFactor(finestResidual);
      }
#ifdef SW_SLAVE
      else
      {
        scalar temp = 0.0;
        init_MVM_Arrays(&arrays1, nCells);
        arrays1.A1Ptr = finestResidual.values();
        arrays1.A2Ptr = b.values();
        arrays1.A3Ptr = Apsi.values();
        arrays1.k1Ptr = &temp;
        residualNormFactor_host(&arrays1);
        reduceSum(&temp);
        solverPerf.finalResidual() = sqrt(temp);
      }
#endif

#ifdef DEBUG
      scalar convergenceRate =
          solverPerf.finalResidual() / solverPerf.previousResidual();
      solverPerf.previousResidual() = solverPerf.finalResidual();
      // swTimer::endTimer("MG Vcycle");
      IFPRINT
      {
        commcator_->log() << "At cycle = ";
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
  }

  return solverPerf;
}
