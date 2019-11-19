#include "unapMatrix.hpp"

//- no preconditioner
void UNAP::matrix::preconditioner::precondition
(
	scalarField       &w,
    const scalarField &r
) const
{
	const label nCells  = w.size();
	const scalar *rPtr  = r.begin();

	forAll(i, nCells)
	{
		w[i] = rPtr[i];
	}
}
