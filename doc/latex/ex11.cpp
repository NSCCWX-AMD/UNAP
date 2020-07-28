#include <math.h>
#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include "readFromOpenFOAM.hpp"
#include <sstream>
#include <string.h>
#include "printUNAP.hpp"

#ifdef SW_SLAVE
#include "swAthread.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

//- test using data output from OpenFOAM

using namespace UNAP;

#define LOCATEFILE(newName, fileName, dir) \
{ \
	std::ostringstream os; \
	os << dir << NPROCS << "/" << fileName << "_" << MYID << ".txt"; \
	strcpy(newName, os.str().c_str()); \
} \

int main()
{
	/* Initialize MPI */
   	unapMPI::initMPI();

   	//- initialize athread environment
#ifdef SW_SLAVE
	swacc_init();
#endif

   	const char* dir = "./exData/openfoam/cavity/20w/p";
   	char fileName[200];

   	if(PARRUN)
   	{
   		UNAP::unapMPI::unapCommunicator().barrier() ;
   	}

   
	UNAPCOUT << "Start reading data" << ENDL;


   	lduMatrix lduA;
   	LOCATEFILE(fileName, "A_p", dir);
   	constructLDUMatrixFromOpenFOAM(lduA, fileName);

   	if(PARRUN)
   	{
   		LOCATEFILE(fileName, "interfaces_p", dir);
   		constructLDUInterfacesFromOpenFOAM(lduA, fileName);
   	}

   	label nCells = lduA.size();
   	scalarField b(nCells);

   	LOCATEFILE(fileName, "b_p", dir);
   	constructVectorFromOpenFOAM(b, fileName);

   	if(PARRUN)
   	{
   		UNAP::unapMPI::unapCommunicator().barrier() ;
   	}

 
	UNAPCOUT << "Finish reading data" << ENDL;
	

   	scalar tol = 0.0;
	scalar relTol = 1e-6;
	label  nFaces = lduA.upper().size();


	const bool useMG = true;
	const bool usePBiCGStab = false;
	scalarField x(nCells, 0.0);

	if(useMG)
	{
		scalarField weights(nFaces);
		forAll(i, nFaces)
		{
			weights[i] = mag(lduA.upper()[i]);
		}

		//- construct coarse grid using upper coefficients
		lduAgglomeration aggl(lduA);
		aggl.agglomerate(weights);
		PtrList<matrix::smoother> sm(aggl.size());

		//- using Gauss-Seidel smoother
		//- be noted that GS is not compatible with MLB
		// forAll(i, aggl.size())
		// {
		// 	lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
		// 	sm.setLevel(i, *smLocPtr);
		// }

		//- using Chebyshev smoother
		forAll(i, aggl.size())
		{
			chebySmoother* smLocPtr = new chebySmoother;
			sm.setLevel(i, *smLocPtr);
		}

		//- construct AMG solver
		MGSolver MG(lduA, aggl, sm);

		//- this part will using MLB to reorder matrix, b and x
#ifdef SW_SLAVE
		lduA.constructMLBIterator();
		lduA.reorderVector(b);
#ifdef SWTIMER
    	swTimer::startTimer("MLB reorder");
#endif
		lduA.reorderVector(x);
#ifdef SWTIMER
    	swTimer::endTimer("MLB reorder");
#endif
		aggl.agglomerationReorderTopo();
		lduA.reorderLDUValues();
#endif

		//- MG controls
		MG.SET_tolerance(tol);
		MG.SET_relTol(relTol);
		MG.SET_nPreSweeps(1);
		MG.SET_maxIter(15);
		MG.SET_ifPrint(true);

#ifdef SWTIMER
    	swTimer::startTimer("MG Solve");
#endif
    	//- solve phase
		matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

#ifdef SW_SLAVE
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

		
		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	else if(usePBiCGStab)
	{
		lduDiagPrecond precond(lduA);

		// lduDICPrecond precond(lduA);

		// lduDILUPrecond precond(lduA);

		PBiCGStab PBiCGStabSolver(precond);

		// PBiCGStabSolver.SET_minIter(5);

		PBiCGStabSolver.SET_maxIter(10);
		PBiCGStabSolver.SET_ifPrint(true);

		matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);


		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	else
	{
		// lduDiagPrecond precond(lduA);

		lduDICPrecond precond(lduA);

		PCG PCGSolver(precond);


		// PCGSolver.SET_minIter(5);

		// PCGSolver.SET_maxIter(50);
		PCGSolver.SET_ifPrint(true);

		matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}

	// printVector(x, "xEnd");

#ifdef SW_SLAVE
	swacc_end();
#endif

#ifdef SWTIMER
	swTimer::printTimer();
#endif

	// MPI_Finalize();
}
