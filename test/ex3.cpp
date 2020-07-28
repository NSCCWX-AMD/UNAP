#include "lduMatrix.hpp"
#include "PCG.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"

//- test for PCG

using namespace UNAP;

int main()
{
	//- A = |3 2|
	//-     |2 6|

	//- b = |2 -8|

	label nCells = 2;
	scalarField x(nCells);

	scalar *bVal = new scalar[2];
	bVal[0] =  2.0;
	bVal[1] = -8.0;
	scalarField b(bVal, nCells);

	labelField lowerAddr(1, 0);
	labelField upperAddr(1, 1);

	scalarField lower(1, 2);
	scalarField upper(1, 2);
	scalarField diag (nCells);

	diag[0] = 3.0;
	diag[1] = 6.0;

	const lduMatrix lduA
	(
		nCells,
	    lowerAddr,
	    upperAddr,
	    lower,
	    diag,
	    upper
	);

	const matrix *A = &lduA;

	// lduDiagPrecond precond(lduA);

	lduDICPrecond precond(lduA);

	// lduDILUPrecond precond(lduA);

	PCG PCGSolver(precond);

	matrix::solverPerformance solverPerf = PCGSolver.solve(x, *A, b);

	UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged to: " << ENDL;
	UNAPCOUT << "x0 = " << x[0] << ", x1 = " << x[1] << ENDL;

	return 0;
}
