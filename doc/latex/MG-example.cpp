//- using upper coefficients in matrix as the weights of coarsening in AMG
//- alternative using face areas
scalarVector weights(nFaces);
forAll(i, nFaces)
{
	weights[i] = mag(lduA.upper()[i]);
}

//- MG setup phase
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
lduA.reorderVector(x);
aggl.agglomerationReorderTopo();
lduA.reorderLDUValues();
#endif

//- MG controls
MG.SET_tolerance(tol);	//- set absolute tolerance
MG.SET_relTol(relTol);	//- set relative tolerance
MG.SET_nPreSweeps(1);   //- pre-smooth numbers in V-cycle
MG.SET_maxIter(15);     //- maximum iteration numbers
MG.SET_ifPrint(true);   //- print information when calculating

#ifdef SWTIMER
swTimer::startTimer("MG Solve");
#endif
//- solve phase
matrix::solverPerformance solverPerf = MG.solve(x, lduA, b);

#ifdef SW_SLAVE
lduA.restoreVector(x);
#endif

#ifdef SWTIMER
swTimer::endTimer("MG Solve");
#endif

//- print iteration numbers

UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;