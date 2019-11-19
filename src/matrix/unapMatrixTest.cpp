#include "unapMatrix.hpp"

bool UNAP::matrix::solverPerformance::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance,
    const label  iter,
    const label  minIter
)
{
    if
    (
     (   finalResidual_ < Tolerance
     || (
            RelTolerance > SMALL
         && finalResidual_ <= RelTolerance*initialResidual_
        )
     )
     &&
     iter >= minIter
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}

bool UNAP::matrix::solverPerformance::checkConvergence
(
    const scalar Tolerance,
    const label  iter,
    const label  minIter
)
{
    if (	(finalResidual_ < Tolerance)
    	&& 	iter >= minIter
       )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


bool UNAP::matrix::solverPerformance::checkConvergence
(
    const scalar Tolerance,
    const scalar RelTolerance
)
{
    if
    (
     (   finalResidual_ < Tolerance
     || (
            RelTolerance > SMALL
         && finalResidual_ <= RelTolerance*initialResidual_
        )
     )
    )
    {
        converged_ = true;
    }
    else
    {
        converged_ = false;
    }

    return converged_;
}


bool UNAP::matrix::solverPerformance::checkSingularity
(
    const scalar residual
)
{
    if (residual > SMALL)
    {
        singular_ = false;
    }
    else
    {
        singular_ = true;
    }

    return singular_;
}

