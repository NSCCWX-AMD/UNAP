#ifndef MG_HPP
#define MG_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class MGSolver
:
	public matrix::solver
{
protected:

    //- finest matrix pointer
    const matrix& finestMatrix_;

    //- if the agglomeration is always stored
	bool cacheAgglomeration_;

	 //- number of pre-smoothing sweeps
    label nPreSweeps_;

    //- number of post-smoothing sweeps
    label nPostSweeps_;

    //- number of smoothing sweeps on finest mesh
    label nFinestSweeps_;

    //- choose if the corrections should be scaled.
    //  by default corrections for symmetric matrices are scaled
    //  but not for asymmetric matrices.
    bool scaleCorrection_;

    //- the agglomeration
    matrix::agglomeration& agglomeration_;

    //- smoother in each level
    // const matrix::smoother& smoother_;
    const PtrList<matrix::smoother>& smoothers_;

    //- calculate and return the scaling factor from Acf, coarseSource
    //  and coarseField.
    //  at the same time do a Jacobi iteration on the coarseField using
    //  the Acf provided after the coarseField values are used for the
    //  scaling factor.
    scalar scalingFactor
    (
        scalarField& Acf,
        scalarField& field,
        const scalarField& source,
        const matrix& A
    ) const;

    //- initialize the data structures for the V-cycle
    void initVcycle
    (
        PtrList<scalarField>& coarseCorrFields,
        PtrList<scalarField>& coarseSources
    ) const;

    //- perform a single GAMG V-cycle with pre, post and finest smoothing.
    void Vcycle
    (
        scalarField& psi,
        const scalarField& source,
        scalarField& Apsi,
        scalarField& finestCorrection,
        scalarField& finestResidual,
        PtrList<scalarField>& coarseCorrFields,
        PtrList<scalarField>& coarseSources
    ) const;

    //- relTol on coarsest level

    //- solve the coarsest level with an iterative solver
    void solveCoarsestLevel
    (
        scalarField& coarsestCorrField,
        const scalarField& coarsestSource
    ) const;

public:

    //- constructors
    MGSolver
    (
        const matrix& A,
        matrix::agglomeration& agglomerator,
        const PtrList<matrix::smoother>& smoothers
    );

	//- destructor
    virtual ~MGSolver();

    //- solve
    virtual matrix::solverPerformance solve
    (
        scalarField& x,
        const matrix& A,
        const scalarField& b
    ) const;

    //- access to parameters
    void SET_nPreSweeps(const label nPreSweeps);

    void SET_nPostSweeps(const label nPostSweeps);

    void SET_nFinestSweeps(const label nFinestSweeps);

    void SET_cacheAgglomeration(bool cacheAgglomeration);

    void SET_scaleCorrection(bool scaleCorrection);

    void initSmoothers();

    matrix::agglomeration& agglomeration()
    {
        return agglomeration_;
    }

};

} //- namespace UNAP

#endif //- MG_HPP
