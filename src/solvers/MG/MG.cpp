#include "MG.hpp"

UNAP::MGSolver::MGSolver
(
	const matrix& A,
	matrix::agglomeration& MGAgglomeration,
	const PtrList<matrix::smoother>& smoothers
)
:
	finestMatrix_(A),
	cacheAgglomeration_(false),
    nPreSweeps_(0),
    nPostSweeps_(2),
    nFinestSweeps_(2),
    scaleCorrection_(A.symm()),
    // scaleCorrection_(false),
    agglomeration_(MGAgglomeration),
    smoothers_(smoothers)
{}

UNAP::MGSolver::~MGSolver()
{
	if(!cacheAgglomeration_)
	{}
}


void UNAP::MGSolver::SET_nPreSweeps(const label nPreSweeps)
{
	nPreSweeps_ = nPreSweeps;
}

void UNAP::MGSolver::SET_nPostSweeps(const label nPostSweeps)
{
	nPostSweeps_ = nPostSweeps;
}

void UNAP::MGSolver::SET_nFinestSweeps(const label nFinestSweeps)
{
	nFinestSweeps_ = nFinestSweeps;
}

void UNAP::MGSolver::SET_cacheAgglomeration(bool cacheAgglomeration)
{
	cacheAgglomeration_ = cacheAgglomeration;
}

void UNAP::MGSolver::SET_scaleCorrection(bool scaleCorrection)
{
	scaleCorrection_ = scaleCorrection;
}


void UNAP::MGSolver::initSmoothers()
{
	const label size = agglomeration_.size();

	forAll(i, size)
	{
		smoothers_[i].init();
	}
}






