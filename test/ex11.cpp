#include <math.h>
#include <string.h>

#include <sstream>

#include "MG.hpp"
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "chebySmoother.hpp"
#include "lduAgglomeration.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include "lduDiagPrecond.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "printUNAP.hpp"
#include "rcmf.h"
#include "readFromOpenFOAM.hpp"

#ifdef SW_SLAVE
#include "swAthread.h"
#include "swlu.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

#ifdef SW_SLAVE
// #define UNAT_MLB
// #define UNAT_RSS
#endif

//- test using data output from OpenFOAM

using namespace UNAP;

#define LOCATEFILE(newName, fileName, dir)                           \
  {                                                                  \
    std::ostringstream os;                                           \
    os << dir << NPROCS << "/" << fileName << "_" << MYID << ".txt"; \
    strcpy(newName, os.str().c_str());                               \
  }

int main()
{
  /* Initialize MPI */
  unapMPI::initMPI();
#ifdef SW_SLAVE
  swacc_init();
#endif

  const char *dir =
      "/home/export/online3/amd_dev1/guhf/UNAP_test_data/exData/openfoam/"
      "cavity/2.5w/p";
  char fileName[200];

  if (PARRUN)
  {
    UNAP::unapMPI::unapCommunicator().barrier();
  }

  UNAPCOUT << "Start reading diagonal of A" << ENDL;

  lduMatrix lduA;
  LOCATEFILE(fileName, "A_p", dir);
  constructLDUMatrixFromOpenFOAM(lduA, fileName);

  UNAPCOUT << "End reading diagonal of A" << ENDL;

  if (PARRUN)
  {
    UNAPCOUT << "Start reading interfaces" << ENDL;

    LOCATEFILE(fileName, "interfaces_p", dir);
    constructLDUInterfacesFromOpenFOAM(lduA, fileName);

    UNAPCOUT << "End reading interfaces" << ENDL;
  }

  label nCells = lduA.size();
  label nFaces = lduA.upper().size();
  scalarVector b(nCells);

  UNAPCOUT << "Start reading b" << ENDL;

  LOCATEFILE(fileName, "b_p", dir);
  constructVectorFromOpenFOAM(b, fileName);

  UNAPCOUT << "End reading b" << ENDL;

  if (PARRUN)
  {
    UNAP::unapMPI::unapCommunicator().barrier();
  }

  UNAPCOUT << "Finish reading data" << ENDL;

  scalar tol = 0.0;
  scalar relTol = 1e-6;

  const bool useMG = false;
  const bool usePBiCGStab = false;
  scalarVector x(nCells, 0.0);

#ifdef SW_SLAVE
  swlu_prof_init();
#endif

  labelVector postV(nCells);
  labelVector postE(nFaces);
  rcmLDU_nowrite(nFaces,
                 nCells,
                 lduA.lowerAddr().begin(),
                 lduA.upperAddr().begin(),
                 postV.begin(),
                 postE.begin());

  if (useMG)
  {
    UNAPCOUT << " *************************************************************"
                " \n\n ";
    UNAPCOUT << "                        use  MG   solver                      "
                " \n\n ";
    UNAPCOUT << " *************************************************************"
                " \n\n ";

    scalarVector weights(nFaces);
    forAll(i, nFaces) { weights[i] = mag(lduA.upper()[i]); }

    lduAgglomeration aggl(lduA);
    aggl.SET_maxLevels(50);
    aggl.agglomerate(weights);
    PtrList<matrix::smoother> sm(aggl.size());

    forAll(i, aggl.size())
    {
      UNAPCOUT << "At coarse level " << i << ":" << ENDL;
      lduMatrix &cm = aggl.coarseMatrix(i);
      label cnCells = cm.size();
      label cnFaces = cm.upperAddr().size();

      UNAPCOUT << "nCells = " << cnCells << ", nFaces = " << cnFaces << ENDL;

      labelVector cpostV(cnCells);
      labelVector cpostE(cnFaces);

      labelVector clowerAddrOld(cm.lowerAddr());
      labelVector cupperAddrOld(cm.upperAddr());

      rcmLDU_nowrite(cnFaces,
                     cnCells,
                     cm.lowerAddr().begin(),
                     cm.upperAddr().begin(),
                     cpostV.begin(),
                     cpostE.begin());
    }

    // forAll(i, aggl.size())
    // {
    // 	lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
    // 	sm.setLevel(i, *smLocPtr);
    // }

    forAll(i, aggl.size())
    {
      chebySmoother *smLocPtr = new chebySmoother;
      sm.setLevel(i, *smLocPtr);
    }
    MGSolver MG(lduA, aggl, sm);

#ifdef UNAT_MLB
    lduA.constructMLBIterator();
    lduA.reorderVector(b);
    lduA.reorderVector(x);
    aggl.agglomerationReorderTopo();
    lduA.reorderLDUValues();
#endif

#ifdef UNAT_RSS
    lduA.constructRSSIterator();
    aggl.agglConstructRSSIterator();
#endif

    MG.SET_tolerance(tol);
    MG.SET_relTol(relTol);
    MG.SET_nPreSweeps(2);
    MG.SET_minIter(1);
    MG.SET_maxIter(1);
    MG.SET_ifPrint(true);

#ifdef SWTIMER
    swTimer::startTimer("MG Solve first cycle");
#endif
    matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);
#ifdef SWTIMER
    swTimer::endTimer("MG Solve first cycle");
#endif

    MG.SET_minIter(30);
    MG.SET_maxIter(30);

#ifdef SWTIMER
    swTimer::startTimer("MG Solve");
#endif

#ifdef SW_SLAVE
    swlu_prof_start();
#endif
    solverPerf = MG.solve(x, lduA, b);
#ifdef SW_SLAVE
    swlu_prof_stop();
#endif

#ifdef UNAT_MLB
#ifdef SWTIMER
    swTimer::startTimer("MLB restore");
#endif
    lduA.restoreVector(x);
#ifdef SWTIMER
    swTimer::endTimer("MLB restore");
#endif
#endif

#ifdef SWTIMER
    swTimer::endTimer("MG Solve");
#endif

    UNAPCOUT << "After " << solverPerf.nIterations()
             << " iterations, the solution is converged!" << ENDL;
    UNAPCOUT << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }
  else if (usePBiCGStab)
  {
    UNAPCOUT << " *************************************************************"
                " \n\n ";
    UNAPCOUT << "                        use  PBiCGStab  solver                "
                " \n\n ";
    UNAPCOUT << " *************************************************************"
                " \n\n ";

#ifdef UNAT_MLB
    lduA.constructMLBIterator();
    lduA.reorderVector(b);
    lduA.reorderVector(x);
    lduA.reorderLDUValues();
#endif

    lduDiagPrecond precond(lduA);

    // lduDICPrecond precond(lduA);

    // lduDILUPrecond precond(lduA);

    PBiCGStab PBiCGStabSolver(precond);

    PBiCGStabSolver.SET_minIter(20);

    PBiCGStabSolver.SET_maxIter(20);
    PBiCGStabSolver.SET_ifPrint(true);

#ifdef UNAT_RSS
    lduA.constructRSSIterator();
#endif

#ifdef SWTIMER
    swTimer::startTimer("PBiCGStab Solve");
#endif

    matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);

#ifdef SWTIMER
    swTimer::endTimer("PBiCGStab Solve");
#endif

#ifdef UNAT_MLB
#ifdef SWTIMER
    swTimer::startTimer("MLB restore");
#endif
    lduA.restoreVector(x);
#ifdef SWTIMER
    swTimer::endTimer("MLB restore");
#endif
#endif

    UNAPCOUT << "After " << solverPerf.nIterations()
             << " iterations, the solution is converged!" << ENDL;
    UNAPCOUT << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }
  else
  {
    UNAPCOUT << " *************************************************************"
                " \n\n ";
    UNAPCOUT << "                          use  PCG  solver                    "
                " \n\n ";
    UNAPCOUT << " *************************************************************"
                " \n\n ";

    // lduDiagPrecond precond(lduA);

    lduDICPrecond precond(lduA);

    PCG PCGSolver(precond);

    // PCGSolver.SET_minIter(5);

    // PCGSolver.SET_maxIter(50);
    PCGSolver.SET_ifPrint(true);

    matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

    UNAPCOUT << "After " << solverPerf.nIterations()
             << " iterations, the solution is converged!" << ENDL;
    UNAPCOUT << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }

  // test scalar byte and label byte

  UNAPCOUT << "test byte: \n";
  UNAPCOUT << "label " << sizeof(label) << " ,scalar " << sizeof(scalar)
           << std::endl;

#ifdef SW_SLAVE
  swlu_prof_print();
#endif
  // printVector(x, "xEnd");

#ifdef SW_SLAVE
  swacc_end();
#endif

#ifdef SWTIMER
  swTimer::printTimer();
#endif
  UNAP::unapMPI::unapCommunicator().barrier();
  // UNAP::unapMPI::exitMPI();
  // MPI_Finalize();
}
