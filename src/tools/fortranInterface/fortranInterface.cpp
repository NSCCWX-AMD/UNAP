#include "unap.hpp"
#include "fortranInterface.hpp"
#include "matrixConversion.hpp"

#include "PCG.hpp"
#include "PBiCGStab.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"

#include "printUNAP.hpp"

// using namespace UNAP;

#define RETURN_VALUE(x, objectX, nCells) \
forAll(i, nCells) \
{ \
	x[i] = objectX[i]; \
}

#define PTR2OBJ(ptr, cls, obj) \
cls* ptr##_tmp = (cls*)*ptr; \
CHECK_POINTER(ptr##_tmp) \
cls& obj = *ptr##_tmp;


void UNAP::ldumatrixcreat_
(
	long int*   APtrPtr,
	label*     nCellsPtr,
	label*     upperSizePtr,
	label*     lowerAddrValue,
	label*     upperAddrValue,
	scalar*    lowerValue,
	scalar*    diagValue,
	scalar*    upperValue
)
{
	label nCells = *nCellsPtr;
	label upperSize = *upperSizePtr;
	labelField lowerAddr(lowerAddrValue, upperSize);
	labelField upperAddr(upperAddrValue, upperSize);

	scalarField diag(diagValue, nCells);
	scalarField upper(upperValue, upperSize);

	scalarField& lower = upper;
	scalarField* lowerPtr = NULL;

	if(lowerValue != upperValue)
	{
		lowerPtr = new scalarField(lowerValue, upperSize);
		lower = *lowerPtr;
	}

	*(lduMatrix**)APtrPtr = new lduMatrix
	(
		nCells,
	    lowerAddr,
	    upperAddr,
	    lower,
	    diag,
	    upper
	);

	DELETE_OBJECT_POINTER(lowerPtr)
}


void UNAP::coo2ldumatrixcreat_
(
	long int*     APtrPtr,
	const scalar* dataPtr,
	const label*  fRowsPtr,
	const label*  fColsPtr,
	const label*  nCellsPtr,
	const label*  sizePtr,
	const label*  symmPtr
)
{
	unapMPI::initMPI();
	label* cRowsPtr = new label[*sizePtr];
	label* cColsPtr = new label[*sizePtr];
	scalar* cValsPtr = new scalar[*sizePtr];

	forAll(i, *sizePtr)
	{
		cRowsPtr[i] = fRowsPtr[i] - 1;
		cColsPtr[i] = fColsPtr[i] - 1;
		cValsPtr[i] = dataPtr[i];
	}

	bool symm;
	if(*symmPtr)
		symm = true;
	else
		symm = false;

	reorderCOO(cValsPtr, cRowsPtr, cColsPtr, *nCellsPtr, *sizePtr);

	lduMatrix& lduA = coo2ldu
	(
		cValsPtr,
		cRowsPtr,
		cColsPtr,
	   *nCellsPtr,
	   *sizePtr,
	    symm
	);

	*(lduMatrix**)APtrPtr = &lduA;

	DELETE_POINTER(cRowsPtr)
	DELETE_POINTER(cColsPtr)
	DELETE_POINTER(cValsPtr)
}


void UNAP::csr2ldumatrixcreat_
(
	long int*     APtrPtr,
	const scalar* dataPtr,
	const label*  fRowsPtr,
	const label*  fColsPtr,
	const label*  nCellsPtr,
	const label*  sizePtr,
	const label*  symmPtr
)
{
	unapMPI::initMPI();
	lduMatrix& lduA = csr2ldu
	(
		dataPtr,
		fRowsPtr,
		fColsPtr,
	   *nCellsPtr,
	   *sizePtr,
	   *symmPtr
	);

	*(lduMatrix**)APtrPtr = &lduA;
}


void UNAP::matrixinterfacescreat_
(
	long int*      APtrPtr,
	const label*   nNeiProcsPtr,
	const label*   destRankPtr,
	const label*   locPositionPtr,
	const label*   faceCellsPtr,
	const scalar*  dataPtr
)
{
	lduMatrix* APtr = (lduMatrix*)*APtrPtr;
	CHECK_POINTER(APtr)
	lduMatrix& lduA = *APtr;

	const label nNeiProcs = *nNeiProcsPtr;

	PtrList<patch>* patchesPtr = new PtrList<patch>(nNeiProcs);

	forAll(intI, nNeiProcs)
	{
		const label neighborID = destRankPtr[intI] - 1;
		const label localSize = locPositionPtr[intI+1] - locPositionPtr[intI];
		patch* patchIPtr = new patch(localSize, MYID, neighborID);

		scalar* localData = new scalar[localSize];
		label*  localFaceCells = new label[localSize];

		forAll(faceI, localSize)
		{
			label start = locPositionPtr[intI] + faceI - 1;
			localData[faceI] = dataPtr[start];
			localFaceCells[faceI] = faceCellsPtr[start] - 1;
		}

		scalarField* patchCoeffsPtr = new scalarField(localData, localSize);
   		labelField* locFaceCellsPtr = new labelField(localFaceCells, localSize);

   		delete [] localData;
   		delete [] localFaceCells;

   		patchIPtr->patchCoeffs(*patchCoeffsPtr);
   		patchIPtr->faceCells(*locFaceCellsPtr);

   		patchesPtr->setLevel(intI, *patchIPtr);
	}

	interfaces* interfacesLocalPtr = new interfaces(*patchesPtr);
	lduA.matrixInterfaces(*interfacesLocalPtr);
}


void UNAP::pcgsolversolve_
(
	scalar*   xValue,
	long int* APtrPtr,
	scalar*   bValue,
	label*    nCellsPtr,
	label*    precondTypePtr,
	scalar*   tolPtr,
	scalar*   relTolPtr,
	label*    maxIterPtr,
	label*    minIterPtr,
	label*    num_iterationsPtr,
	scalar*   res_normPtr
)
{
	label nCells = *nCellsPtr;
	scalarField x(xValue, nCells);
	scalarField b(bValue, nCells);

	lduMatrix* APtr = (lduMatrix*)*APtrPtr;
	CHECK_POINTER(APtr)
	lduMatrix& lduA = *APtr;

	label precondType = *precondTypePtr;

	matrix::preconditioner* precondPtr = NULL;

	if(precondType == 1)
	{
		precondPtr = new lduDiagPrecond(lduA);
	}
	else if(precondType == 2)
	{
		precondPtr = new lduDICPrecond(lduA);
	}
	else if(precondType == 3)
	{
		precondPtr = new lduDILUPrecond(lduA);
	}

	PCG PCGSolver(*precondPtr);
	PCGSolver.SET_maxIter(*maxIterPtr);
	PCGSolver.SET_minIter(*minIterPtr);
	PCGSolver.SET_tolerance(*tolPtr);
	PCGSolver.SET_relTol(*relTolPtr);

	matrix::solverPerformance solverPerf = PCGSolver.solve(x, lduA, b);

#ifdef DEBUG

	if(solverPerf.converged())
	{
		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	else
	{
		UNAPCOUT << "The PCG solver reaches the maximum iterations." << ENDL;
	}

#endif

	RETURN_VALUE(xValue, x, nCells)
	*num_iterationsPtr = solverPerf.nIterations();
	*res_normPtr = solverPerf.finalResidual();

	DELETE_OBJECT_POINTER(precondPtr)
	DELETE_OBJECT_POINTER(APtr)
}


void UNAP::pbicgstabsolversolve_
(
	scalar*   xValue,
	long int* APtrPtr,
	scalar*   bValue,
	label*    nCellsPtr,
	label*    precondTypePtr,
	scalar*   tolPtr,
	scalar*   relTolPtr,
	label*    maxIterPtr,
	label*    minIterPtr,
	label*    num_iterationsPtr,
	scalar*   res_normPtr
)
{
	label nCells = *nCellsPtr;
	scalarField x(xValue, nCells);
	scalarField b(bValue, nCells);

	lduMatrix* APtr = (lduMatrix*)*APtrPtr;
	CHECK_POINTER(APtr)
	lduMatrix& lduA = *APtr;

	label precondType = *precondTypePtr;

	matrix::preconditioner* precondPtr = NULL;

	if(precondType == 1)
	{
		precondPtr = new lduDiagPrecond(lduA);
	}
	else if(precondType == 2)
	{
		precondPtr = new lduDICPrecond(lduA);
	}
	else if(precondType == 3)
	{
		precondPtr = new lduDILUPrecond(lduA);
	}

	PBiCGStab PBiCGStabSolver(*precondPtr);

	PBiCGStabSolver.SET_maxIter(*maxIterPtr);
	PBiCGStabSolver.SET_minIter(*minIterPtr);
	PBiCGStabSolver.SET_tolerance(*tolPtr);
	PBiCGStabSolver.SET_relTol(*relTolPtr);
	PBiCGStabSolver.SET_ifPrint(true);

	// printLDUMatrix(lduA, "old_A_u");
	// printInterfaces(lduA, "old_interfaces_u");
	// printVector(b, "old_b");
	// return;

	matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);

#ifdef DEBUG

	if(solverPerf.converged())
	{
		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	else
	{
		UNAPCOUT << "The PBiCGStab solver reaches the maximum iterations." << ENDL;
	}

#endif

	RETURN_VALUE(xValue, x, nCells)
	*num_iterationsPtr = solverPerf.nIterations();
	*res_normPtr = solverPerf.finalResidual();

	DELETE_OBJECT_POINTER(precondPtr)
	DELETE_OBJECT_POINTER(APtr)
}


void UNAP::mgsolversolve_
(
	scalar*   xValue,
	long int* APtrPtr,
	scalar*   bValue,
	label*    nCellsPtr,
	label*    agglTypePtr,
	label*    smootherTypePtr,
	scalar*   tolPtr,
	scalar*   relTolPtr,
	label*    maxIterPtr,
	label*    minIterPtr,
	label*    num_iterationsPtr,
	scalar*   res_normPtr,
	scalar*   faceAreaPtr
)
{
	label nCells = *nCellsPtr;
	scalarField x(xValue, nCells);
	scalarField b(bValue, nCells);

	lduMatrix* APtr = (lduMatrix*)*APtrPtr;
	CHECK_POINTER(APtr)
	lduMatrix& lduA = *APtr;

	// printLDUMatrix(lduA, "unap_A_p");
	// printInterfaces(lduA, "unap_interfaces_p");

	UNAP::unapMPI::unapCommunicator().barrier();

	label agglType = *agglTypePtr;

	matrix::agglomeration* agglPtr = NULL;

	if(agglType == 1)
	{
		scalarField weights(lduA.upper().size());

		scalar average = lduA.upper().SumSqrt();

		forAll(i, weights.size())
		{
			weights[i] = mag(lduA.upper()[i])/average;
		}
		agglPtr = new lduAgglomeration(lduA);
		agglPtr->agglomerate(weights);
	}
	else if(agglType == 2)
	{
		scalarField weights(faceAreaPtr, lduA.upper().size());
		agglPtr = new lduAgglomeration(lduA);

		// printVector(weights, "old_facearea");
		// return;

		agglPtr->agglomerate(weights);
	}

	label smootherType = *smootherTypePtr;
	PtrList<matrix::smoother> sm(agglPtr->size());

	if(smootherType == 1)
	{
		forAll(i, (*agglPtr).size())
		{
			lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
			sm.setLevel(i, *smLocPtr);
		}
	}
	else if(smootherType == 2)
	{
		forAll(i, (*agglPtr).size())
		{
			chebySmoother* smLocPtr = new chebySmoother;
			sm.setLevel(i, *smLocPtr);
		}
	}

	MGSolver MG(lduA, *agglPtr, sm);

	MG.SET_maxIter(*maxIterPtr);
	MG.SET_minIter(*minIterPtr);
	MG.SET_tolerance(*tolPtr);
	MG.SET_relTol(*relTolPtr);

	MG.SET_nPreSweeps(0);
	MG.SET_ifPrint(true);

	// printLDUMatrix(lduA, "old_A_p");
	// printInterfaces(lduA, "old_interfaces_p");
	// printVector(b, "old_b");
	// return;

	matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

#ifdef DEBUG

	if(solverPerf.converged())
	{
		UNAPCOUT << "After " << solverPerf.nIterations() << " cycles, the solution is converged!" << ENDL;
	}
	else
	{
		UNAPCOUT << "The Multigrid solver reaches the maximum cycles." << ENDL;
	}

#endif

	RETURN_VALUE(xValue, x, nCells)
	*num_iterationsPtr = solverPerf.nIterations();
	*res_normPtr = solverPerf.finalResidual();

	DELETE_OBJECT_POINTER(agglPtr)
	DELETE_OBJECT_POINTER(APtr)
}


void UNAP::reordercoo_
(
	scalar* val,
	label*  row,
	label*  col,
	label*  nCellsPtr,
	label*  sizePtr
)
{
	if(row[0] != 0)
	{
		forAll(i, *sizePtr)
		{
			row[i] = row[i] - 1;
			col[i] = col[i] - 1;
		}
	}
	reorderCOO(val, row, col, *nCellsPtr, *sizePtr);
}


void UNAP::reorderuface__
(
	label*  row,
	label*  col,
	label*  nCellsPtr,
	label*  sizePtr,
	label*  newOrder
)
{
	forAll(i, *sizePtr)
	{
		row[i] = row[i] - 1;
		col[i] = col[i] - 1;
	}

	reorderUFace(row, col, *nCellsPtr, *sizePtr, newOrder);
}


void UNAP::reorderlface__
(
	label*  row,
	label*  col,
	label*  nCellsPtr,
	label*  sizePtr,
	label*  newOrder
)
{
	forAll(i, *sizePtr)
	{
		row[i] = row[i] - 1;
		col[i] = col[i] - 1;
	}

	reorderLFace(row, col, *nCellsPtr, *sizePtr, newOrder);
}


void UNAP::reordervalue__
(
	scalar* val,
	label*  newOrder,
	label*  sizePtr
)
{
	reorderValue(val, newOrder, *sizePtr);
}


void UNAP::contruct_sw_matrix__
(
	long int*      APtrPtr,
	const label*   nCellsPtr,
	const label*   rowsPtr,
	const label*   colsPtr,
	const label*   sizePtr
)
{
	unapMPI::initMPI();

	const label nCells = *nCellsPtr;
	const label size = *sizePtr;

	labelField* upperAddrPtr = new labelField(size);
	labelField& upperAddr = *upperAddrPtr;

	labelField* lowerAddrPtr = new labelField(size);
	labelField& lowerAddr = *lowerAddrPtr;

	forAll(i, size)
	{
		upperAddr[i] = colsPtr[i];
		lowerAddr[i] = rowsPtr[i];
	}

	scalarField* diagPtr = new scalarField(nCells);
	scalarField& diag = *diagPtr;

	scalarField* upperPtr = new scalarField(size);
	scalarField& upper = *upperPtr;

	lduMatrix* lduAPtr = new lduMatrix(nCells, lowerAddr, upperAddr, upper, diag, upper, true);
	*(lduMatrix**)APtrPtr = lduAPtr;
}


void UNAP::contruct_sw_matrix_interfaces__
(
	long int*     APtrPtr,
	const label*  nNeiProcsPtr,
	const label*  destRankPtr,
	const label*  offDiagRowsPtr,
	const label*  offDiagStartsPtr
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	const label nNeiProcs = *nNeiProcsPtr;

	PtrList<patch>* patchesPtr = new PtrList<patch>(nNeiProcs);

	forAll(intI, nNeiProcs)
	{
		const label neighborID = destRankPtr[intI] - 1;
		const label localSize = offDiagStartsPtr[intI+1] - offDiagStartsPtr[intI];
		patch* patchIPtr = new patch(localSize, MYID, neighborID);

		label* localFaceCells = new label[localSize];

		forAll(faceI, localSize)
		{
			label start = offDiagStartsPtr[intI] + faceI - 1;
			localFaceCells[faceI] = offDiagRowsPtr[start] - 1;
		}

   		labelField* locFaceCellsPtr = new labelField(localFaceCells, localSize, true);
   		scalarField* localDataPtr = new scalarField(localSize);

   		patchIPtr->faceCells(*locFaceCellsPtr);
   		patchIPtr->patchCoeffs(*localDataPtr);

   		patchesPtr->setLevel(intI, *patchIPtr);
	}

	interfaces* interfacesLocalPtr = new interfaces(*patchesPtr);
	lduA.matrixInterfaces(*interfacesLocalPtr);
}

#ifdef SW_SLAVE
void UNAP::construct_mlb_iterator__
(
	long int*     APtrPtr
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)
	lduA.constructMLBIterator();
}
#endif


void UNAP::fill_sw_matrix_coefficients__
(
	long int*      APtrPtr,
	const scalar*  diagPtr,
	const scalar*  upperPtr,
	const scalar*  lowerPtr
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	bool symm = upperPtr == lowerPtr ? true : false;

	const label nCells = lduA.size();
	const label nFaces = lduA.upperAddr().size();

	// UNAPCOUT << "In unap, nCells = " << nCells << ", nFaces = " << nFaces << ENDL;

	scalar* upperData = lduA.upper().begin();
	scalar* diagData = lduA.diag().begin();

	forAll(i, nCells)
	{
		diagData[i] = diagPtr[i];
	}

// #ifdef SW_SLAVE
// 	lduA.reorderVector(lduA.diag());
// #endif

	forAll(i, nFaces)
	{
		upperData[i] = upperPtr[i];
	}

	if(!symm)
	{
		// UNAPCOUT << "symm = " << symm << ENDL;
		scalarField lower(lowerPtr, nFaces);
		lduA.SET_lower(lower);
	}
	else if(!lduA.symm())
	{
		lduA.freeLower();
	}

#ifdef SW_SLAVE
	lduA.reorderLDUValues();
#endif
}


void UNAP::fill_sw_matrix_interfaces_coefficients__
(
	long int*      APtrPtr,
	const label*   offDiagStartsPtr,
	const scalar*  offDiagCoeffs
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	interfaces& lduAInter = lduA.matrixInterfaces();

	const label nNeiProcs = lduAInter.size();

	forAll(intI, nNeiProcs)
	{
		patch& patchI = lduAInter.patchList(intI);
		const label localSize = offDiagStartsPtr[intI+1] - offDiagStartsPtr[intI];
		scalar* localData = patchI.patchCoeffs().begin();

		forAll(faceI, localSize)
		{
			label start = offDiagStartsPtr[intI] + faceI - 1;
			localData[faceI] = offDiagCoeffs[start];
		}
	}
}


//- MultiGrid solver solve and controls
void UNAP::contruct_solver_mg__
(
	long int*      mgPtrPtr,
	long int*      APtrPtr,
	long int*      AgglPtr,
	const scalar*  weightsPtr,
	const label*   smootherTypePtr,
	const label*   maxLevelsPtr,
	const label*   nCellsCoarsestPtr
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	const label nFaces = lduA.upper().size();

	scalarField weights(weightsPtr, nFaces);

	// printVector(weights, "new_facearea");
	// UNAPCOUT << "finish writing" << ENDL;
	// return;

	lduAgglomeration* agglPtr = new lduAgglomeration(lduA);

	agglPtr->SET_maxLevels(*maxLevelsPtr);
	agglPtr->SET_nCellsInCoarsestLevel(*nCellsCoarsestPtr);
	agglPtr->agglomerate(weights);

	const label coarseLevels = agglPtr->size();

	PtrList<matrix::smoother>* smoothersPtr = new PtrList<matrix::smoother> (coarseLevels);

	const label smootherType = *smootherTypePtr;
	PtrList<matrix::smoother>& sm = *smoothersPtr;

	if(smootherType == 1)
	{
		
		UNAPCOUT << "Gauss-Seidel smoother used in AMG." << ENDL;

		forAll(i, coarseLevels)
		{
			lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
			sm.setLevel(i, *smLocPtr);
		}
	}
	else if(smootherType == 2)
	{
		
		UNAPCOUT << "Chebyshev smoother used in AMG." << ENDL;
		
		forAll(i, coarseLevels)
		{
			chebySmoother* smLocPtr = new chebySmoother;
			sm.setLevel(i, *smLocPtr);
		}
	}
	else
	{
		UNAPCOUT << "ERROR: smoother type is not recognized!" << ENDL;
		UNAPCOUT << "Valid option is: 1--GaussSeidel, 2--Chebyshev" << ENDL;
	}

	MGSolver* MGPtr = new MGSolver(lduA, *agglPtr, sm);

	*(MGSolver**)mgPtrPtr = MGPtr;

	*(lduAgglomeration**)AgglPtr = agglPtr;
}

#ifdef SW_SLAVE
void UNAP::mg_coarse_mlb__(long int* agglPtr)
{
	PTR2OBJ(agglPtr, lduAgglomeration, aggl);
	aggl.agglomerationReorderTopo();
}
#endif


void UNAP::sw_solver_mg_set_maxiter__
(
	long int*    solverPtrPtr,
	const label* maxIterPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver);
	solver.SET_maxIter(*maxIterPtr);
}


void UNAP::sw_solver_mg_set_miniter__
(
	long int*    solverPtrPtr,
	const label* minIterPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver)
	solver.SET_minIter(*minIterPtr);
}


void UNAP::sw_solver_mg_set_tol__
(
	long int*     solverPtrPtr,
	const scalar* tolPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver)
	solver.SET_tolerance(*tolPtr);
}


void UNAP::sw_solver_mg_set_reltol__
(
	long int*     solverPtrPtr,
	const scalar* reltolPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver)
	solver.SET_relTol(*reltolPtr);
}


void UNAP::sw_solver_mg_set_npresweeps__
(
	long int*     solverPtrPtr,
	const label*  numPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver);
	solver.SET_nPreSweeps(*numPtr);
}


void UNAP::sw_solver_mg_set_npostsweeps__
(
	long int*     solverPtrPtr,
	const label*  numPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver)
	solver.SET_nPostSweeps(*numPtr);
}


void UNAP::sw_solver_mg_set_nfinestsweeps__
(
	long int*     solverPtrPtr,
	const label*  numPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver);
	solver.SET_nFinestSweeps(*numPtr);
}


void UNAP::sw_solve_mg__
(
	long int*      mgPtrPtr,
	long int*      APtrPtr,
	scalar*        xPtr,
	scalar*        bPtr,
	scalar* 	   res_normPtr
)
{
	PTR2OBJ(mgPtrPtr, MGSolver, MG)
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	MG.initSmoothers();

	const label nCells = lduA.size();

	scalarField x(xPtr, nCells);
	scalarField b(bPtr, nCells);

	MG.SET_ifPrint(true);

	// printLDUMatrix(lduA, "new_A_p");
	// printInterfaces(lduA, "new_interfaces_p");
	// printVector(b, "new_b");
	// UNAPCOUT << "finish writing" << ENDL;
	// return;

#ifdef SW_SLAVE
	lduA.reorderVector(b);
	lduA.reorderVector(x);
#endif

	printMessage("Start solving in AMG solver");
	matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

#ifdef SW_SLAVE
	lduA.restoreVector(x);
#endif

	RETURN_VALUE(xPtr, x, nCells);

#ifdef DEBUG
	
	
	if(solverPerf.converged())
	{
		UNAPCOUT << "After " << solverPerf.nIterations() << " cycles, the solution is converged!" << ENDL;
	}
	else
	{
		UNAPCOUT << "The Multigrid solver reaches the maximum cycles." << ENDL;
	}
	
#endif

	*res_normPtr = solverPerf.finalResidual();
}


//- PBiCGStab solver solve and controls
void UNAP::contruct_solver_pbicgstab__
(
	long int*      solverPtrPtr
)
{
	PBiCGStab* solverPtr = new PBiCGStab();

	*(PBiCGStab**)solverPtrPtr = solverPtr;
}

void UNAP::sw_solver_pbicgstab_set_maxiter__
(
	long int*    solverPtrPtr,
	const label* maxIterPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver);
	solver.SET_maxIter(*maxIterPtr);
}


void UNAP::sw_solver_pbicgstab_set_miniter__
(
	long int*    solverPtrPtr,
	const label* minIterPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver)
	solver.SET_minIter(*minIterPtr);
}


void UNAP::sw_solver_pbicgstab_set_tol__
(
	long int*     solverPtrPtr,
	const scalar* tolPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver);
	solver.SET_tolerance(*tolPtr);
}


void UNAP::sw_solver_pbicgstab_set_reltol__
(
	long int*     solverPtrPtr,
	const scalar* reltolPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver);
	solver.SET_relTol(*reltolPtr);
}


void UNAP::sw_solver_pbicgstab_set_precond__
(
	long int*     solverPtrPtr,
	long int*     APtrPtr,
	const label*  precondTypePtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver)
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	const label type = *precondTypePtr;
	if(type == 1)
	{
		lduDiagPrecond* pptr1 = new lduDiagPrecond(lduA);
		solver.SET_preconditioner(*pptr1);
		printMessage("Diagonal preconditioner is used");
	}
	else if(type == 2)
	{
		lduDICPrecond* pptr1 = new lduDICPrecond(lduA);
		solver.SET_preconditioner(*pptr1);
		printMessage("DIC preconditioner is used");
	}
	else if(type == 3)
	{
		lduDILUPrecond* pptr1 = new lduDILUPrecond(lduA);
		solver.SET_preconditioner(*pptr1);
		printMessage("DILU preconditioner is used");
	}
	else
	{
		UNAPCOUT << "ERROR: preconditioner type is not recognized!" << ENDL;
		UNAPCOUT << "Valid option is: 1--diagonal, 2--DIC, 3--DILU" << ENDL;
	}
}


void UNAP::sw_solve_pbicgstab__
(
	long int*      solverPtrPtr,
	long int*      APtrPtr,
	scalar*        xPtr,
	scalar*        bPtr,
	scalar* 	   res_normPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver)
	PTR2OBJ(APtrPtr, lduMatrix, lduA)

	const label nCells = lduA.size();

	scalarField x(xPtr, nCells);
	scalarField b(bPtr, nCells);

	solver.SET_ifPrint(true);

	

#ifdef SW_SLAVE
	lduA.reorderVector(b);
	lduA.reorderVector(x);
#endif

	printMessage("Start solving in PBiCGStab solver");
	matrix::solverPerformance solverPerf = solver.solve(x, lduA, b);

#ifdef SW_SLAVE
	lduA.restoreVector(x);
#endif

	RETURN_VALUE(xPtr, x, nCells);

#ifdef DEBUG
	
	if(solverPerf.converged())
	{
		UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	else
	{
		UNAPCOUT << "The PBiCGStab solver reaches the maximum iterations." << ENDL;
	}
	
#endif

	*res_normPtr = solverPerf.finalResidual();
}


void UNAP::sw_matrix_destroy__
(
	long int* APtrPtr
)
{
	PTR2OBJ(APtrPtr, lduMatrix, lduA)
	delete &lduA;
}


void UNAP::sw_solver_destroy_mg__
(
	long int* solverPtrPtr
)
{
	PTR2OBJ(solverPtrPtr, MGSolver, solver)
	delete &solver;
}


void UNAP::sw_solver_destroy_pbicgstab__
(
	long int* solverPtrPtr
)
{
	PTR2OBJ(solverPtrPtr, PBiCGStab, solver)
	delete &solver;
}

