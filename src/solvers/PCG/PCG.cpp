#include "PCG.hpp"

#define IFPRINT \
if(!MYID && ifPrint_)

UNAP::PCG::PCG()
:
	deletePrecondPtr_(false),
	precondPtr_(NULL)
{
	precondPtr_ = new matrix::preconditioner;
	deletePrecondPtr_ = true;
}

UNAP::PCG::PCG
(
	matrix::preconditioner& precond
)
:
	deletePrecondPtr_(false),
	precondPtr_(NULL)
{
	if(&precond == NULL)
	{
		precondPtr_ = new matrix::preconditioner;
		deletePrecondPtr_ = true;
	}
	else
	{
		precondPtr_ = &precond;
	}
}


UNAP::matrix::solverPerformance UNAP::PCG::solve
(
	scalarField& x,
	const matrix& A,
	const scalarField& b
) const
{
	matrix::solverPerformance solverPerf;

	label nCells = x.size();
	scalar* xPtr = x.values();
	const scalar* bPtr = b.values();

	scalarField pA(nCells);
	scalar* pAPtr = pA.values();

	scalarField wA(nCells);
	scalar* wAPtr = wA.values();

	scalarField rA(nCells);
	scalar* rAPtr = rA.values();

	scalar wArA = GREAT;
	scalar wArAold = wArA;

	//- calculate A.psi
	A.spMV(wA, x);

	//- calculate initial residual field
	forAll(i, nCells)
	{
		rAPtr[i] = bPtr[i] - wAPtr[i];
	}

#ifdef DEBUG
	//- calculate normalisation factor
	scalar normFactor = this->normFactor(b);
#endif

	solverPerf.initialResidual() = this->normFactor(rA);
	solverPerf.finalResidual() = solverPerf.initialResidual();

#ifdef DEBUG
IFPRINT
{
	COUT << "At nIter = ";
	std::cout.width(5);
	COUT << solverPerf.nIterations();
	COUT << ",   ini res = ";
	std::cout.width(11);
	COUT << solverPerf.initialResidual();
	COUT << ",   rel res = ";
	std::cout.width(11);
	COUT << solverPerf.initialResidual()/solverPerf.initialResidual();
	COUT << ",   b norm = ";
	std::cout.width(11);
	std::cout.setf(std::ios::scientific);
	COUT << normFactor << ENDL;
}
#endif


	do
	{
		//- store previous wArA
		wArAold = wArA;

		precondPtr_->precondition(wA, rA);

		//- update search directions
		wArA = dot(wA, rA);

        // --- Test for singularity
        if (solverPerf.checkSingularity(mag(wArA)))
        {
#ifdef DEBUG
            IFPRINT
            {
            	COUT << "singularity! wArA = " << wArA << ENDL;
            }
#endif
            break;
        }

		if(solverPerf.nIterations() == 0)
		{
			forAll(i, nCells)
			{
				pAPtr[i] = wAPtr[i];
			}
		}
		else
		{
			scalar beta = wArA/wArAold;

			forAll(i, nCells)
			{
				pAPtr[i] = wAPtr[i] + beta*pAPtr[i];
			}
		}

		//- update preconditioned residual
		A.spMV(wA, pA);

		scalar wApA = dot(wA, pA);

		//- update solution and residual
		scalar alpha = wArA/wApA;

		forAll(i, nCells)
		{
			xPtr[i]  += alpha*pAPtr[i];
			rAPtr[i] -= alpha*wAPtr[i];
		}

		// solverPerf.finalResidual() = rA.SumMag() / normFactor;
		solverPerf.finalResidual() = this->normFactor(rA);

#ifdef DEBUG
IFPRINT
{
		COUT << "At nIter = ";
		std::cout.width(5);
		COUT << solverPerf.nIterations()+1;
		COUT << ",   fin res = ";
		std::cout.width(11);
		COUT << solverPerf.finalResidual();
		COUT << ",   rel res = ";
		std::cout.width(11);
		COUT << solverPerf.finalResidual()/normFactor << ENDL;
}
#endif
	} while
	(
		(++solverPerf.nIterations() < maxIter_
    &&
        !(solverPerf.checkConvergence
         	(
              	tolerance_,
              	relTol_,
              	solverPerf.nIterations(),
              	minIter_
         	)
 		 )
        )
	);

	return solverPerf;
}
