#include <math.h>
#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "PBiCGStab.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"

//-test for matrix conversion: coo2ldu
//- using MultiGrid solver to solve a 2D Poisson equation

//- problem:
//- (d^2)u/d(x^2) + (d^2)u/d(y^2) = -8.0*pi*pi*sin(2*pi*x)*sin(2*pi*y)
//- analytical solution is: sin(2*pi*x)*sin(2*pi*y)

#define N  10
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
	label myid, neighborid, num_procs;
   	MPI_Init(&argc, &argv);

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

	// forAll(i, size)
	// {
	// 	COUT << "At i = " << i << ", data is " << data[i] << ENDL;
	// }

	bool symm = true;

	lduMatrix& lduA = coo2ldu
	(
		data,
		row,
		col,
		nCells,
		size,
		symm
	);

	// lduAgglomeration aggl(lduA, lduA.upper());
	// PtrList<matrix::smoother> sm(aggl.size());

	// forAll(i, aggl.size())
	// {
	// 	lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
	// 	sm.setLevel(i, *smLocPtr);
	// }

	// forAll(i, aggl.size())
	// {
	// 	chebySmoother* smLocPtr = new chebySmoother;
	// 	sm.setLevel(i, *smLocPtr);
	// }
	// MGSolver MG(lduA, aggl, sm);

	// MG.SET_tolerance(1e-10);
	// MG.SET_nPreSweeps(4);

	// matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

	PBiCGStab PBiCGStabSolver;
	PBiCGStabSolver.SET_maxIter(1);
	PBiCGStabSolver.SET_minIter(1);

	matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);

	COUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;

	delete []data;
	delete []col;
	delete []row;

	delete &lduA;

	return 0;
}
