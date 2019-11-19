#include "chebySmoother.hpp"
#include "eigenDiagPCG.hpp"
#include "lduDiagPrecond.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

void UNAP::chebySmoother::smooth
(
	scalarField       &x,
    const lduMatrix   &A,
    const scalarField &b,
    const label       nSweeps
) const
{
	lduDiagPrecond precond(A);
	if(eigFirstTimeComputed_)
	{
		eigenDiagPCG calEigen
		(
			A,
			x,
			b,
			precond,
			nDiagPCGs_
		);
		maxEigPCG_ = calEigen.maxEigenValue();
		eigFirstTimeComputed_ = false;

#ifdef DEBUG
		// if(!MYID)
		// {
		// 	std::cout.setf(std::ios::scientific);
		// 	COUT << "Factor is " << maxEigPCG_ << ENDL;
		// }
#endif
	}
	else
	{
#ifdef SWTIMER
    	swTimer::startTimer("cheby pure computing");
#endif

		register const label nCells = x.size();
	    scalar* xPtr = x.begin();

	    scalarField p(nCells);
	    scalar* pPtr = p.begin();

	    scalarField w(nCells);
	    scalar* wPtr = w.begin();

    	const scalar* rDPtr = precond.rD().begin();

	    const scalar* bPtr = b.begin();

		scalar maxEig, minEig;
	    maxEig = maxEigPCG_ * boostFactorCheby_;
	    minEig = maxEigPCG_ / eigRatioCheby_;

	    scalar delta, theta, s1;
	    scalar rhok, rhokp1, dtemp1, dtemp2;
	    delta = 2.0 / (maxEig - minEig);
	    theta = 0.5 * (maxEig + minEig);
	    s1    = theta * delta;
	    rhok  = 1.0 / s1;

	    //- calculate A.x
	    A.spMV(w, x);

#ifdef SW_SLAVE
		MVM_Arrays arrays1;
#endif
	     // At sweep = 0, update x
        scalar Rtheta = 1.0 / theta;
        IFNOT_SWACC
        {
            forAll(cell, nCells)
            {
                wPtr[cell]  = rDPtr[cell] * (bPtr[cell] - wPtr[cell]) * Rtheta;
                xPtr[cell] += wPtr[cell];
            }
        }
#ifdef SW_SLAVE
        else
        {
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr    = wPtr;
            arrays1.A2Ptr    = xPtr;
            arrays1.A3Ptr    = (scalar*)rDPtr;
            arrays1.A4Ptr    = (scalar*)bPtr;
            arrays1.k1       = Rtheta;
            arrays1.returnA2 = 1;
            // wPtr  = rDPtr * (bPtr - wPtr) * Rtheta
            // xPtr += wPtr
            vectorOps_host(&arrays1, &slave_userFunc_aEcMuSdMiaSMuk1_bEbPa);
        }
#endif

	    for (label sweep=1; sweep<nSweeps; sweep++)
	    {
	        rhokp1 = 1.0 / (2.0 * s1 - rhok);
	        dtemp1 = rhokp1 * rhok;
	        dtemp2 = 2.0 * rhokp1 * delta;
	        rhok   = rhokp1;

	        // Timer::startTimer("Cheby Amul");
	        A.spMV(p, x);
	        // Timer::endTimer("Cheby Amul");

	        IFNOT_SWACC
            {
                forAll(cell, nCells)
                {
                    wPtr[cell]  = dtemp1 * wPtr[cell] + dtemp2 * rDPtr[cell] * (bPtr[cell] - pPtr[cell]);
                    xPtr[cell] += wPtr[cell];
                }
            }
#ifdef SW_SLAVE
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = wPtr;
                arrays1.A2Ptr = (scalar*)rDPtr;
                arrays1.A3Ptr = (scalar*)bPtr;
                arrays1.A4Ptr = pPtr;
                arrays1.k1    = dtemp1;
                arrays1.k2    = dtemp2;
                // wPtr = dtemp1 * wPtr + dtemp2 * rDPtr * (bPtr - pPtr)
                vectorOps_host(&arrays1, &slave_userFunc_aEk1MuaPk2MubMuScMidS);

                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = xPtr;
                arrays1.A2Ptr = wPtr;
                // xPtr += wPtr
                vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
            }
#endif
	    }
#ifdef SWTIMER
    	swTimer::endTimer("cheby pure computing");
#endif
	}
}


void UNAP::chebySmoother::smooth
(
	scalarField       &x,
    const matrix 	  &A,
    const scalarField &b,
    const label       nSweeps
) const
{
	smooth
	(
		x,
		static_cast<const lduMatrix &> (A),
		b,
		nSweeps
	);
}
