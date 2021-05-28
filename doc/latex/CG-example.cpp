//- construct PBiCGStab solver by matrix A
lduDiagPrecond precond(lduA);

//- preconditioners
//- DIC
// lduDICPrecond precond(lduA);

//- DILU
// lduDILUPrecond precond(lduA);

//- Diagonal(Jacobi)
PBiCGStab PBiCGStabSolver(precond);

//- solver controls
PBiCGStabSolver.SET_minIter(1);    //- set minimum iteration numbers
PBiCGStabSolver.SET_maxIter(10);   //- set maximum iteration numbers
PBiCGStabSolver.SET_ifPrint(true); //- print information when calculating

//- solve phase
matrix::solverPerformance solverPerf = PBiCGStabSolver.solve(x, lduA, b);


UNAPCOUT << "After " << solverPerf.nIterations() << " iterations, the solution is converged!" << ENDL;