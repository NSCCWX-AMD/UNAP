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

#include "rcmf.h"

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
#ifdef SW_SLAVE
	swacc_init();
#endif



   	const char* dir = "./../../../exData/openfoam/cavity/2.5w/p";
   	char fileName[200];

   	if(PARRUN)
   	{
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

   	if(!MYID)
	{
		COUT << "Start reading diagonal of A" << ENDL;
	}

   	lduMatrix lduA;
   	LOCATEFILE(fileName, "A_p", dir);
   	constructLDUMatrixFromOpenFOAM(lduA, fileName);

   	if(!MYID)
	{
		COUT << "End reading diagonal of A" << ENDL;
	}


   	if(PARRUN)
   	{
   		if(!MYID)
		{
			COUT << "Start reading interfaces" << ENDL;
		}
   		LOCATEFILE(fileName, "interfaces_p", dir);
   		constructLDUInterfacesFromOpenFOAM(lduA, fileName);
   		if(!MYID)
		{
			COUT << "End reading interfaces" << ENDL;
		}
   	}

   	label nCells = lduA.size();
   	label nFaces = lduA.upper().size();
   	scalarField b(nCells);

   	if(!MYID)
	{
		COUT << "Start reading b" << ENDL;
	}
   	LOCATEFILE(fileName, "b_p", dir);
   	constructVectorFromOpenFOAM(b, fileName);
   	if(!MYID)
	{
		COUT << "End reading b" << ENDL;
	}

   	if(PARRUN)
   	{
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

   	if(!MYID)
	{
		COUT << "Finish reading data" << ENDL;
	}

   	scalar tol = 0.0;
	scalar relTol = 1e-6;



	const bool useMG = true;
	const bool usePBiCGStab = false;
	scalarField x(nCells, 0.0);

#ifdef SW_SLAVE
	swlu_prof_init();
#endif

	labelField postV(nCells);
   	labelField postE(nFaces);
   	rcmLDU_nowrite(nFaces, nCells, lduA.lowerAddr().begin(), lduA.upperAddr().begin(), postV.begin(), postE.begin());


	if(useMG)
	{
		if(!MYID){
			COUT <<" ************************************************************* \n\n ";
			COUT <<"                        use  MG   solver                       \n\n ";
			COUT <<" ************************************************************* \n\n ";
		}
		scalarField weights(nFaces);
		forAll(i, nFaces)
		{
			weights[i] = mag(lduA.upper()[i]);
		}

		lduAgglomeration aggl(lduA);
		aggl.SET_maxLevels(50);
		aggl.agglomerate(weights);
		PtrList<matrix::smoother> sm(aggl.size());

		forAll(i, aggl.size())
		{
			COUT << "At coarse level " << i << ":" << ENDL;
			lduMatrix& cm = aggl.coarseMatrix(i);
			label cnCells = cm.size();
		   	label cnFaces = cm.upperAddr().size();

		   	COUT << "nCells = " << cnCells << ", nFaces = " << cnFaces << ENDL;

		   	labelField cpostV(cnCells);
		   	labelField cpostE(cnFaces);

		   	labelField clowerAddrOld(cm.lowerAddr());
		   	labelField cupperAddrOld(cm.upperAddr());

		   	rcmLDU_nowrite(cnFaces, cnCells, cm.lowerAddr().begin(), cm.upperAddr().begin(), cpostV.begin(), cpostE.begin());
		}


		// forAll(i, aggl.size())
		// {
		// 	lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
		// 	sm.setLevel(i, *smLocPtr);
		// }

		forAll(i, aggl.size())
		{
			chebySmoother* smLocPtr = new chebySmoother;
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

		if(!MYID)
		{
			COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
			COUT << "finalResidual " << solverPerf.finalResidual()  << ENDL;
		}

	}
	else if(usePBiCGStab)
	{
		if(!MYID){
			COUT <<" ************************************************************* \n\n ";
			COUT <<"                        use  PBiCGStab  solver                 \n\n ";
			COUT <<" ************************************************************* \n\n ";
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

		if(!MYID)
		{
			COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
			COUT << "finalResidual " << solverPerf.finalResidual()  << ENDL;
		}
	}
	else
	{
		if(!MYID){
			COUT <<" ************************************************************* \n\n ";
			COUT <<"                          use  PCG  solver                     \n\n ";
			COUT <<" ************************************************************* \n\n ";
		}
		// lduDiagPrecond precond(lduA);

		lduDICPrecond precond(lduA);

		PCG PCGSolver(precond);


		// PCGSolver.SET_minIter(5);

		// PCGSolver.SET_maxIter(50);
		PCGSolver.SET_ifPrint(true);

		matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

		if(!MYID)
		{
			COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
			COUT << "finalResidual " << solverPerf.finalResidual()  << ENDL;
		}
	}

	// test scalar byte and label byte
	if(!MYID){
		std::cout<<"test byte: \n";
		std::cout<<"label "<<sizeof(label)<<" ,scalar "<<sizeof(scalar)<<std::endl;
	}

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

	MPI_Finalize();
}
