#include "MGPreconditioner.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"

UNAP::MGPrecond::MGPrecond
(
	const lduMatrix& A
)
:
	nVcycles_(2),
	APtr_(&A)
{
	const scalarVector& upper = A.upper();
	scalarVector weights(upper.size());
	forAll(i, weights.size())
	{
		weights[i] = mag(upper[i]);
	}

	lduAgglomeration aggl(A);
	aggl.agglomerate(weights);
	PtrList<matrix::smoother> sm(aggl.size());

	forAll(i, aggl.size())
	{
		lduGaussSeidelSmoother* smLocPtr = new lduGaussSeidelSmoother;
		sm.setLevel(i, *smLocPtr);
	}

	MGSolver MG(A, aggl, sm);
}


void UNAP::MGPrecond::precondition
(
	scalarVector& w,
	const scalarVector& r
) const
{

}
