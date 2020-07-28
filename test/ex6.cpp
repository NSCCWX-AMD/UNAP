#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"

#include "Timers.hpp"

//- test for lduMatrix
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

int main()
{
	label nCells = pow(N-1, 2);
	scalarVector x(nCells, 0.0);

	//- domain solved
	const scalar x0 = 0.0;
	const scalar x1 = 1.0;
	const scalar y0 = 0.0;
	const scalar y1 = 1.0;

	const scalar dx = (x1 - x0) / N;
	const scalar dy = (y1 - y0) / N;

	scalar *bAll = new scalar[(N+1)*(N+1)];
	forAll(i, N+1)
		forAll(j, N+1)
		{
			bAll[i*(N+1)+j] = -8.0*pi*pi*sin(2*pi*i*dx)*sin(2*pi*j*dy);
		}

	// UNAPCOUT << "bAll: " << ENDL;
	// forAll(i, (N+1)*(N+1))
	// {
	// 	UNAPCOUT << bAll[i] << ENDL;
	// }

	scalar *xAll = new scalar[(N+1)*(N+1)];
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

	// UNAPCOUT << "xAll: " << ENDL;
	// forAll(i, (N+1)*(N+1))
	// {
	// 	UNAPCOUT << xAll[i] << ENDL;
	// }

	scalarVector b(nCells);

	scalar rhx2 = 1.0/dx/dx;
	scalar rhy2 = 1.0/dy/dy;

	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			b[(i-1)*(N-1) + j -1] = bAll[i*(N+1)+j] - Au(xAll, rhx2, rhy2, i, j);
		}

	label nZeros = 0;
	labelVector nZerosCells(nCells+1);
	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			label n0 = nZeros;
			if((i+1)<N) nZeros++;
			if((j+1)<N) nZeros++;
			label n1 = nZeros;

			nZerosCells[(i-1)*(N-1) + j - 1 + 1] = n1 - n0 + nZerosCells[(i-1)*(N-1) + j - 1];
		}

	labelVector lowerAddr(nZeros);
	labelVector upperAddr(nZeros);

	UNAPCOUT << "Number of non zeros in upper coefficients is " << nZeros << ENDL;

	label nStart = 0;
	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			if((j+1)<N)
			{
				lowerAddr[nStart] = (i-1)*(N-1) + j - 1;
				upperAddr[nStart] = (i-1)*(N-1) + j - 1 + 1;
				nStart++;
			}

			if((i+1)<N)
			{
				lowerAddr[nStart] = (i-1)*(N-1) + j - 1;
				upperAddr[nStart] = (i-1)*(N-1) + j - 1 + N - 1;
				nStart++;
			}
		}

	scalarVector upper(nZeros, rhx2);
	scalarVector &lower = upper;
	scalarVector diag (nCells,   -4*rhx2);

	// exit(0);

	const lduMatrix lduA
	(
		nCells,
	    lowerAddr,
	    upperAddr,
	    lower,
	    diag,
	    upper
	);

	scalar tol = 1e-10;
	scalar relTol = 1e-6;

	scalarVector weights(upper.size());
	forAll(i, weights.size())
	{
		weights[i] = mag(upper[i]);
	}

	lduAgglomeration aggl(lduA);
	aggl.agglomerate(upper);
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

	MG.SET_tolerance(tol);
	MG.SET_relTol(relTol);
	MG.SET_nPreSweeps(4);

	matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

	UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;
	UNAPCOUT << "Let me check now: " << ENDL;

	label okNums = 0;
	scalar relErr = 0.001;
	for(int i=1; i<N; i++)
		for(int j=1; j<N; j++)
		{
			scalar xPos = j * dx;
			scalar yPos = i * dy;
			label pos = (i-1)*(N-1) + j - 1;

			scalar uExact = sin(2*pi*xPos)*sin(2*pi*yPos);
			scalar err = fabs((x[pos] - uExact) / (uExact + SMALL));

			if(err > relErr && fabs(uExact) > max(relTol, tol))
				UNAPCOUT << "err = " << err
					 << ", uExact = " << uExact
					 << ", uComput = " << x[pos] << ENDL;
			else
				okNums++;
		}

	if(okNums == nCells)
	{
		UNAPCOUT << "No cell's error is larger than " << relErr*100 << "%." << ENDL;
	}
	else
	{
		UNAPCOUT << "The number of cells is " << nCells
			 << ", while only " << okNums
			 << " cells have correct solutions." << ENDL;
	}

	delete [] bAll;
	delete [] xAll;

	Timer::printTimer();

	return 0;
}
