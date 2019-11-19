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

//-test for mpi
//- using MultiGrid solver to solve a 2D Poisson equation

//- problem:
//- (d^2)u/d(x^2) + (d^2)u/d(y^2) = -8.0*pi*pi*sin(2*pi*x)*sin(2*pi*y)
//- analytical solution is: sin(2*pi*x)*sin(2*pi*y)

#define N  1000
#define pi 3.1415926535897931

using namespace UNAP;

scalar Au(scalar *u, scalar rhx2, scalar rhy2, label i, label j)
{
	return(
			rhx2*(u[(i-1)*(N+1) + j  ] - 2*u[i*(N+1) + j] + u[(i+1)*(N+1) + j  ]) +
			rhy2*(u[    i*(N+1) + j-1] - 2*u[i*(N+1) + j] + u[    i*(N+1) + j+1])
		  );
}

int main(int argc, char *argv[])
{
	/* Initialize MPI */
   	MPI_Init(&argc, &argv);
   	// unapMPI::init(argc, argv);
   	unapMPI::init();
   	label myid = MYID;
   	label num_procs = NPROCS;

	//- domain solved
	const scalar x0 = 0.0;
	const scalar x1 = 1.0;
	const scalar y0 = 0.0;
	const scalar y1 = 1.0;

	const scalar dx = (x1 - x0) / N;
	const scalar dy = (y1 - y0) / N;

	const scalar rhx2 = 1.0/dx/dx;
	const scalar rhy2 = 1.0/dy/dy;

	label nCells = pow(N-1, 2);

	scalar*  bAll = new scalar[(N+1)*(N+1)];
	forAll(i, N+1)
		forAll(j, N+1)
		{
			bAll[i*(N+1)+j] = -8.0*pi*pi*sin(2*pi*i*dx)*sin(2*pi*j*dy);
		}

	scalar* xAll = new scalar[(N+1)*(N+1)];
	forAll(i, N+1)
		forAll(j, N+1)
		{
			if(	   (i==0)
				|| (i==N)
				|| (j==0)
				|| (j==N))
			{
				xAll[i*(N+1)+j] = sin(2*pi*i*dx)*sin(2*pi*j*dy);
			}
			else
			{
				xAll[i*(N+1)+j] = 0.0;
			}
		}

	scalarField b(nCells);
	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			b[(i-1)*(N-1) + j -1] = bAll[i*(N+1)+j] - Au(xAll, rhx2, rhy2, i, j);
		}

	delete []bAll;
	delete []xAll;

	scalarField x(nCells, 0.0);

	label nZeros = 0;
	labelField nZerosCells(nCells+1);
	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			label n0 = nZeros;
			if((i+1)<N) nZeros++;
			if((j+1)<N) nZeros++;
			label n1 = nZeros;

			nZerosCells[(i-1)*(N-1) + j - 1 + 1] = n1 - n0 + nZerosCells[(i-1)*(N-1) + j - 1];
		}

	label size = 2 * nZeros + nCells;
	label ndx = N - 1;
	label size0 = 0;

	label*  row  = new label[size];
	label*  col  = new label[size];
	scalar* data = new scalar[size];

	for(int i=0; i<ndx; i++)
		for(int j=0; j<ndx; j++)
		{
			int position = i * ndx + j;

			//- down
			if(i > 0)
			{
				int down = position - ndx;
				data[size0] = rhx2;
				row[size0] = position;
				col[size0] = down;
				size0++;
			}

			//- left
			if(j > 0)
			{
				int left = position - 1;
				data[size0] = rhx2;
				row[size0] = position;
				col[size0] = left;
				size0++;
			}

			//- position
			data[size0] = -4*rhx2;
			row[size0] = position;
			col[size0] = position;
			size0++;

			//- right
			if(j < ndx -1)
			{
				int right = position + 1;
				data[size0] = rhx2;
				row[size0] = position;
				col[size0] = right;
				size0++;
			}

			//- up
			if(i < ndx -1)
			{
				int up = position + ndx;
				data[size0] = rhx2;
				row[size0] = position;
				col[size0] = up;
				size0++;
			}
		}

	COUT << "global nnz = " << size << ENDL;

   	label local_size, extra;
   	label ilower, iupper;

   	local_size = nCells / num_procs;
   	extra = nCells - local_size*num_procs;

   	ilower = local_size*myid;
   	ilower += min(myid, extra);

   	iupper = local_size*(myid+1);
   	iupper += min(myid+1, extra);
   	iupper = iupper - 1;

   	/* How many rows do I have? */
   	local_size = iupper - ilower + 1;

   	scalarField local_b(local_size);
   	scalarField local_x(local_size);

   	forAll(i, local_size)
   	{
   		local_b[i] = b[i+ilower];
   		local_x[i] = x[i+ilower];
   	}

   	label local_nnz = 0;

   	for(int i=0; i<size; i++)
   	{
   		label rowI = row[i];
   		if(rowI >= ilower && rowI <= iupper)
   		{
   			local_nnz++;
   		}
   	}

   	COUT << "local_nnz = " << local_nnz << ENDL;

   	scalar* local_data = new scalar[local_nnz];
   	label* local_row = new label[local_nnz];
   	label* local_col = new label[local_nnz];


   	label ith = 0;
   	for(int i=0; i<size; i++)
   	{
   		label rowI = row[i];
   		if(rowI >= ilower && rowI <= iupper)
   		{
   			local_data[ith] = data[i];
   			local_row[ith] = rowI - ilower;
   			local_col[ith] = col[i];
   			ith++;
   		}
   	}

   	label local_diag_size = 0;
   	label local_offdiag_size = 0;
   	for(int i=0; i<local_nnz; i++)
   	{
   		// label rowI = local_row[i];
   		label colI = local_col[i];
   		if(colI >= ilower && colI <=iupper)
   		{
   			local_diag_size++;
   		}
   		else
   		{
   			local_offdiag_size++;
   		}
   	}

   	scalar* local_diag_data = new scalar[local_diag_size];
   	label* local_diag_row = new label[local_diag_size];
   	label* local_diag_col = new label[local_diag_size];

   	scalar* local_offdiag_data = new scalar[local_offdiag_size];
   	label* local_offdiag_row = new label[local_offdiag_size];
   	label* local_offdiag_col = new label[local_offdiag_size];


   	label ndiag = 0;
   	label noffdiag = 0;
   	for(int i=0; i<local_nnz; i++)
   	{
   		// label rowI = local_row[i];
   		label colI = local_col[i];
   		if(colI >= ilower && colI <=iupper)
   		{
   			local_diag_data[ndiag] = local_data[i];
   			local_diag_row[ndiag] = local_row[i];
   			local_diag_col[ndiag] = local_col[i] - ilower;
   			ndiag++;
   		}
   		else
   		{
   			local_offdiag_data[noffdiag] = local_data[i];

   			if(!myid)
   			{
   				local_offdiag_row[noffdiag] = local_row[i];
   				local_offdiag_col[noffdiag] = local_col[i];
   			}

   			noffdiag++;
   		}
   	}

   	label neighborid = myid == 0? 1 : 0;

   	if(!myid && PARRUN)
   	{
   		MPI_Send
	    (
	        &local_offdiag_col[0],
	        local_offdiag_size,
	        MPI_LABEL,
	        neighborid,
	        1,
	        MPI_COMM_WORLD
	    );
   	}

   	if(myid)
   	{
   		MPI_Recv
	    (
	        &local_offdiag_row[0],
	        local_offdiag_size,
	        MPI_LABEL,
	        neighborid,
	        1,
	        MPI_COMM_WORLD,
	        MPI_STATUSES_IGNORE
	    );
   	}

   	patch* patchIPtr = new patch(local_offdiag_size, myid, neighborid);

   	forAll(i, local_offdiag_size)
   	{
   		local_offdiag_row[i] -= ilower;
   	}

   	scalarField* patchCoeffsPtr = new scalarField(local_offdiag_data, local_offdiag_size);
   	labelField* faceCellsPtr = new labelField(local_offdiag_row, local_offdiag_size);

   	patchIPtr->patchCoeffs(*patchCoeffsPtr);
   	patchIPtr->faceCells(*faceCellsPtr);

   	label patchSize = 1;

   	PtrList<patch> patches(patchSize);

   	patches.setLevel(0, *patchIPtr);

   	interfaces* interfacesLocalPtr = new interfaces(patches);

   	//- especially for 2 mpis

	bool symm = true;

	lduMatrix& lduA = coo2ldu
	(
		local_diag_data,
		local_diag_row,
		local_diag_col,
		local_size,
		local_diag_size,
		symm
	);

	lduA.matrixInterfaces(*interfacesLocalPtr);


	//- MG Solver
	lduAgglomeration aggl(lduA);
	aggl.agglomerate(lduA.upper());
	PtrList<matrix::smoother> sm(aggl.size());

	forAll(i, aggl.size())
	{
		lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
		sm.setLevel(i, *smLocPtr);
	}

	// forAll(i, aggl.size())
	// {
	// 	chebySmoother* smLocPtr = new chebySmoother;
	// 	sm.setLevel(i, *smLocPtr);
	// }

	MGSolver MG(lduA, aggl, sm);
	MG.SET_tolerance(1e-10);
	MG.SET_nPreSweeps(4);
	matrix::solverPerformance solverPerf = MG.solve(local_x, lduA, local_b);

	//- PBiCGStab Solver
	// PBiCGStab PBiCGStabSolver;
	// PBiCGStabSolver.SET_maxIter(50);
	// PBiCGStabSolver.SET_minIter(5);
	// matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(local_x, lduA, local_b);

	if(!MYID)
	{
		COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	}
	delete []data;
	delete []col;
	delete []row;

	delete &lduA;

	return 0;
}
