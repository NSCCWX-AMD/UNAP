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
    os << dir << nProcs << "/" << fileName << "_" << myId << ".txt"; \
    strcpy(newName, os.str().c_str());                               \
  }

int main()
{
  /* Initialize MPI */
  unapMPI init;
#ifdef SW_SLAVE
  swacc_init();
#endif
  Communicator *other_comm = init.unapCommunicator();
  label nProcs = other_comm->getMySize();
  label myId = other_comm->getMyId();
  const char *dir = "./../../../exData/openfoam/cavity/2.5w/p";
  char fileName[200];

  if (nProcs > 1)
  {
    other_comm->barrier();
  }

  if (!myId) std::cout << "Start reading diagonal of A" << ENDL;

  lduMatrix lduA(other_comm);
  LOCATEFILE(fileName, "A_p", dir);
  constructLDUMatrixFromOpenFOAM(lduA, fileName);

  if (!myId) std::cout << "End reading diagonal of A" << ENDL;

  if (nProcs > 1)
  {
    if (!myId) std::cout << "Start reading interfaces" << ENDL;

    LOCATEFILE(fileName, "interfaces_p", dir);
    constructLDUInterfacesFromOpenFOAM(lduA, fileName);

    if (!myId) std::cout << "End reading interfaces" << ENDL;
  }

  label nCells = lduA.size();
  label nFaces = lduA.upper().size();
  scalarVector b(nCells, other_comm);

  if (!myId) std::cout << "Start reading b" << ENDL;

  LOCATEFILE(fileName, "b_p", dir);
  constructVectorFromOpenFOAM(b, fileName);

  if (!myId) std::cout << "End reading b" << ENDL;

  if (nProcs > 1)
  {
    other_comm->barrier();
  }

  if (!myId) std::cout << "Finish reading data" << ENDL;

  scalar tol = 0.0;
  scalar relTol = 1e-6;

  const bool useMG = true;
  const bool usePBiCGStab = false;
  scalarVector x(nCells, 0.0, other_comm);

#ifdef SW_SLAVE
  swlu_prof_init();
#endif

  labelVector postV(nCells, other_comm);
  labelVector postE(nFaces, other_comm);
  rcmLDU_nowrite(nFaces,
                 nCells,
                 lduA.lowerAddr().begin(),
                 lduA.upperAddr().begin(),
                 postV.begin(),
                 postE.begin());

  if (useMG)
  {
    if (!myId)
    {
      std::cout
          << " *************************************************************"
             " \n\n ";
      std::cout
          << "                        use  MG   solver                      "
             " \n\n ";
      std::cout
          << " *************************************************************"
             " \n\n ";
    }

    scalarVector weights(nFaces, other_comm);
    forAll(i, nFaces) { weights[i] = mag(lduA.upper()[i]); }

    lduAgglomeration aggl(lduA);
    aggl.SET_maxLevels(50);
    aggl.agglomerate(weights);

    PtrList<matrix::smoother> sm(aggl.size());
    if (myId) std::cout.rdbuf(nullptr);
    forAll(i, aggl.size())
    {
      std::cout << "At coarse level " << i << ":" << ENDL;
      lduMatrix &cm = aggl.coarseMatrix(i);
      label cnCells = cm.size();
      label cnFaces = cm.upperAddr().size();

      std::cout << "nCells = " << cnCells << ", nFaces = " << cnFaces << ENDL;

      labelVector cpostV(cnCells, other_comm);
      labelVector cpostE(cnFaces, other_comm);

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
      chebySmoother *smLocPtr = new chebySmoother(other_comm);
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
    if (!myId)
      std::cout << "After " << solverPerf.nIterations()
                << " iterations, the solution is converged!" << ENDL
                << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }
  else if (usePBiCGStab)
  {
    if (!myId)
    {
      std::cout
          << " *************************************************************"
             " \n\n ";
      std::cout
          << "                        use  PBiCGStab  solver                "
             " \n\n ";
      std::cout
          << " *************************************************************"
             " \n\n ";
    }

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
    if (!myId)
      std::cout << "After " << solverPerf.nIterations()
                << " iterations, the solution is converged!" << ENDL
                << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }
  else
  {
    if (!myId)
    {
      std::cout
          << " *************************************************************"
             " \n\n ";
      std::cout
          << "                          use  PCG  solver                    "
             " \n\n ";
      std::cout
          << " *************************************************************"
             " \n\n ";
    }

    // lduDiagPrecond precond(lduA);

    lduDICPrecond precond(lduA);

    PCG PCGSolver(precond);

    // PCGSolver.SET_minIter(5);

    // PCGSolver.SET_maxIter(50);
    PCGSolver.SET_ifPrint(true);

    matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

    if (!myId)
      std::cout << "After " << solverPerf.nIterations()
                << " iterations, the solution is converged!" << ENDL
                << "finalResidual " << solverPerf.finalResidual() << ENDL;
  }

  // test scalar byte and label byte

  if (!myId)
    std::cout << "test byte: \n"
              << "label " << sizeof(label) << " ,scalar " << sizeof(scalar)
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
  other_comm->barrier();
  // UNAP::unapMPI::exitMPI();
  // MPI_Finalize();
}
