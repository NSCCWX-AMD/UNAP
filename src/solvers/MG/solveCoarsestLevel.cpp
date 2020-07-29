#include "MG.hpp"
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"

//- need to be improved
//- decide preconditioner here

void UNAP::MGSolver::solveCoarsestLevel(
    scalarVector &coarsestCorrField, const scalarVector &coarsestSource) const
{
#ifdef DEBUG1

  UNAPCOUT << "Solving in the coarsest level" << ENDL;

#endif
  const label coarsestLevel = agglomeration_.size() - 1;

  coarsestCorrField = 0;

  matrix::solverPerformance coarseSolverPerf;

  const matrix &coarsestMatrix =
      agglomeration_.coarseMatrixLevels(coarsestLevel);

  scalar relTol = 1.0e-2;
  scalar tol = 0.0;

  if (coarsestMatrix.symm())
  {
    lduDICPrecond precond((lduMatrix &)(coarsestMatrix));
    PCG PCGSolver(precond);

    PCGSolver.SET_tolerance(tol);
    PCGSolver.SET_relTol(relTol);
    PCGSolver.SET_maxIter(100);
    PCGSolver.SET_ifPrint(false);
    coarseSolverPerf =
        PCGSolver.solve(coarsestCorrField, coarsestMatrix, coarsestSource);
  }
  else
  {
    lduDILUPrecond precond((lduMatrix &)(coarsestMatrix));
    PBiCGStab PBiCGStabSolver(precond);

    PBiCGStabSolver.SET_tolerance(tol);
    PBiCGStabSolver.SET_relTol(relTol);
    PBiCGStabSolver.SET_maxIter(100);
    PBiCGStabSolver.SET_ifPrint(false);
    coarseSolverPerf = PBiCGStabSolver.solve(
        coarsestCorrField, coarsestMatrix, coarsestSource);
  }
}
