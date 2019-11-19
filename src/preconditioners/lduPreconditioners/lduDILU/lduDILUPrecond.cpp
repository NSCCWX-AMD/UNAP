#include "lduDILUPrecond.hpp"
#include "lduMatrix.hpp"

UNAP::lduDILUPrecond::lduDILUPrecond
(
	const lduMatrix &A
)
:
	rD_(A.diag())
{
	APtr_ = &A;
	calcReciprocalD(rD_, A);
}

void UNAP::lduDILUPrecond::calcReciprocalD
(
	scalarField     &rD,
	const lduMatrix &A
)
{
	scalar *rDPtr = rD.begin();

    const label  *uPtr      = A.upperAddr().begin();
    const label  *lPtr      = A.lowerAddr().begin();
    const scalar *upperPtr  = A.upper().begin();
    const scalar *lowerPtr  = A.lower().begin();

    // Calculate the DIC diagonal
    const label nFaces = A.upper().size();
    forAll(i, nFaces)
    {
        rDPtr[uPtr[i]] -= upperPtr[i]*lowerPtr[i]/rDPtr[lPtr[i]];
    }


    // Calculate the reciprocal of the preconditioned diagonal
    const label nCells = A.nCells();

    forAll(i, nCells)
    {
        rDPtr[i] = 1.0/rDPtr[i];
    }
}


void UNAP::lduDILUPrecond::precondition
(
	scalarField       &w,
	const scalarField &r
) const
{
    scalar       *wPtr  = w.begin();
    const scalar *rPtr  = r.begin();
    const scalar *rDPtr = rD_.begin();

    const label  *uPtr      = APtr_->upperAddr().begin();
    const label  *lPtr      = APtr_->lowerAddr().begin();
    const scalar *upperPtr  = APtr_->upper().begin();
    const scalar *lowerPtr  = APtr_->lower().begin();

    //- losortAddr has not been defined
    const label  *losortPtr = APtr_->losortAddr().begin();

    label nCells   = w.size();
    label nFaces   = APtr_->upper().size();
    label nFacesM1 = nFaces - 1;

    forAll(i, nCells)
    {
        wPtr[i] = rDPtr[i]*rPtr[i];
    }

    label sface;

    forAll(i, nFaces)
    {
        sface = losortPtr[i];
        wPtr[uPtr[sface]] -=
            rDPtr[uPtr[sface]]*lowerPtr[sface]*wPtr[lPtr[sface]];
    }

    for (label face=nFacesM1; face>=0; face--)
    {
        wPtr[lPtr[face]] -= rDPtr[lPtr[face]]*upperPtr[face]*wPtr[uPtr[face]];
    }
}
