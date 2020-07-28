#include "unapMatrix.hpp"

UNAP::matrix::solver::solver()
:
	maxIter_(1000),
	minIter_(0),
	relTol_(1e-6),
	tolerance_(0.0),
	ifPrint_(false)
{}

scalar UNAP::matrix::solver::normFactor(const scalarField &source) const
{
	// return source.SumMag() + SMALL;
	return source.SumSqrt() + SMALL;
}

void UNAP::matrix::solver::SET_maxIter(label maxIter)
{
    maxIter_ = maxIter;
}

void UNAP::matrix::solver::SET_minIter(label minIter)
{
    minIter_ = minIter;
}

void UNAP::matrix::solver::SET_relTol(scalar relTol)
{
    relTol_ = relTol;
}

void UNAP::matrix::solver::SET_tolerance(scalar tolerance)
{
    tolerance_ = tolerance;
}

void UNAP::matrix::solver::SET_ifPrint(bool ifPrint)
{
    ifPrint_ = ifPrint;
}
