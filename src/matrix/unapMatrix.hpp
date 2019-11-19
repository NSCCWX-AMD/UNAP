#ifndef UNAPMATRIX_HPP
#define UNAPMATRIX_HPP

#include "unapVector.hpp"
#include "interfaces.hpp"

namespace UNAP
{
class matrix
{

public:
	//- destructor
	virtual ~matrix()
	{}

	//- matrix multiplication with updated interfaces
	virtual void spMV
	(
		scalarField       &Apsi,
		const scalarField &psi
	) const = 0;

    //- number of equations( = number of cells to be solved)
    virtual label size() const = 0;

    //- access to diagonal data
    virtual scalarField& diag() const = 0;

    //- matrix is symmetric or asymmetric
    virtual bool symm() const = 0;

    //- access to interfaces
    virtual interfaces& matrixInterfaces() const = 0;

    //- initialize interfaces
    virtual void initInterfaces(const scalarField& psi) const = 0;

    //- update interfaces
    virtual void updateInterfaces(scalarField& Apsi) const = 0;

	//- Class returned by the solver, containing performance statistics
	class solverPerformance
	{
	private:
		scalar initialResidual_;
		scalar finalResidual_;
        scalar previousResidual_;
        label  noIterations_;
        bool   converged_;
        bool   singular_;

    public:
    	//- Construct null
        solverPerformance()
        :
            initialResidual_(0),
            finalResidual_(0),
            previousResidual_(0),
            noIterations_(0),
            converged_(false),
            singular_(false)
        {}

        //- Construct from components
        solverPerformance
        (
            const scalar iRes,
            const scalar fRes,
            const label  nIter,
            const bool   converged,
            const bool   singular
        )
        :
            initialResidual_(iRes),
            finalResidual_(fRes),
            previousResidual_(iRes),
            noIterations_(nIter),
            converged_(converged),
            singular_(singular)
        {}

        // Member functions

        //- Return initial residual
        scalar initialResidual() const
        {
            return initialResidual_;
        }

        //- return initial residual
        scalar& initialResidual()
        {
            return initialResidual_;
        }


        //- return final residual
        scalar finalResidual() const
        {
            return finalResidual_;
        }

        //- return final residual
        scalar& finalResidual()
        {
            return finalResidual_;
        }


        //- return previous residual
        scalar previousResidual() const
        {
            return previousResidual_;
        }

        //- return previous residual
        scalar& previousResidual()
        {
            return previousResidual_;
        }


        //- return number of iterations
        label nIterations() const
        {
            return noIterations_;
        }

        //- return number of iterations
        label &nIterations()
        {
            return noIterations_;
        }


        //- has the solver converged?
        bool converged() const
        {
            return converged_;
        }

        //- is the matrix singular?
        bool singular() const
        {
            return singular_;
        }

        //- convergence test
        bool checkConvergence
        (
            const scalar tolerance,
            const scalar relTolerance,
            const label  iter,
            const label  minIter
        );

        bool checkConvergence
        (
            const scalar tolerance,
            const label  iter,
            const label  minIter
        );

        bool checkConvergence
        (
            const scalar tolerance,
            const scalar relTolerance
        );

        //- singularity test
        bool checkSingularity(const scalar residual);

	};


    //- abstract base-class for matrix agglomeration
    class agglomeration
    {
    protected:

    public:

        //- constructors
        agglomeration();

        //- destructor
        virtual ~agglomeration()
        {}

        //- number of coarse matrix levels
        virtual label size() const = 0;

        //- access to coarse matrix level
        virtual const matrix&  coarseMatrixLevels(const label leveli) const = 0;

        //- access to restrictAddressing_
        //- return cell restrict addressing of given level
        virtual const labelField& restrictAddressing(const label leveli) const = 0;

        //- access to faceRestrictAddressing_
        //- return face restrict addressing of given level
        virtual const labelField& faceRestrictAddressing(const label leveli) const = 0;

        //- agglomerate coarse matrix
        virtual void agglomerateMatrix(const label fineLevelIndex) = 0;

        virtual void agglomerate(const scalarField& weights) = 0;

        virtual void SET_nCellsInCoarsestLevel(const label i) = 0;

        virtual void SET_maxLevels(const label i) = 0;

        // //- temporary
        // virtual void agglomerationReorderTopo() = 0;

        //- restrict (integrate by summation) cell field
        template<typename T>
        void restrictField
        (
            Vector<T>       &cf,
            const Vector<T> &ff,
            const label     fineLevelIndex
        ) const;

        //- prolong (interpolate by injection) cell field
        template<typename T>
        void prolongField
        (
            Vector<T>       &ff,
            const Vector<T> &cf,
            const label     coarseLevelIndex
        ) const;

        //- restrict (integrate by summation) face field
        template<typename T>
        void restrictFaceField
        (
            Vector<T>       &cf,
            const Vector<T> &ff,
            const label     fineLevelIndex
        ) const;
    };


	//- abstract base-class for matrix solvers
	class solver
	{
	protected:

        //- maximum number of iterations in the iteration
        label  maxIter_;

        //- minimum number of iterations in the iteration
        label  minIter_;

        //- convergence tolerance relative to the initial
        scalar relTol_;

        //- final convergence tolerance
        scalar tolerance_;

        //- print converging process
        bool ifPrint_;

    public:

    	//- constructors
    	solver();

    	//- destructor
        virtual ~solver()
        {}

    	virtual solverPerformance solve
		(
			scalarField       &x,
			const matrix 	  &A,
			const scalarField &b
		) const = 0;

		//- return the matrix norm used to normalize the residual for the
        //- stopping criterion
        scalar normFactor
        (
            const scalarField  &source
        ) const;

        inline label maxIter() const
        {
            return maxIter_;
        }

        inline label minIter() const
        {
            return minIter_;
        }

        inline scalar relTol() const
        {
            return relTol_;
        }

        inline scalar tolerance() const
        {
            return tolerance_;
        }

        inline bool ifPrint() const
        {
            return ifPrint_;
        }

        void SET_maxIter(label maxIter);

        void SET_minIter(label minIter);

        void SET_relTol(scalar relTol);

        void SET_tolerance(scalar tolerance);

        void SET_ifPrint(bool ifPrint);
	};

	//- abstract base-class for matrix preconditioners
    class preconditioner
    {
    private:
        scalarField* rDTempPtr_;

    public:
    	//- destructor
        virtual ~preconditioner()
        {}

    	//- return wA the preconditioned form of residual rA
        virtual void precondition
        (
            scalarField       &wA,
            const scalarField &rA
        ) const;

        virtual const scalarField& rD() const
        {
            return *rDTempPtr_;
        }
    };


    //- abstract base-class for matrix smoothers
    class smoother
    {

    public:

        //- constructor
        smoother()
        {}

        //- destructor
        virtual ~smoother()
        {}

        virtual void smooth
        (
            scalarField       &x,
            const matrix      &A,
            const scalarField &b,
            const label       nSweeps
        ) const = 0;

        virtual void init() const = 0;
    };

};


template<typename T>
void matrix::agglomeration::restrictField
(
    Vector<T>       &cf,
    const Vector<T> &ff,
    const label     fineLevelIndex
) const
{
    const labelField &fineToCoarse = restrictAddressing(fineLevelIndex);

#ifdef DEBUG
    if (ff.size() != fineToCoarse.size())
    {

        COUT << "Error in restrictField: field does not correspond to level " << fineLevelIndex
             << " sizes: field = " << ff.size()
             << " level = " << fineToCoarse.size()
             << ENDL;
        ERROR_EXIT;
    }
#endif

    cf.SET_zero();
    label len = ff.size();

    forAll(i, len)
    {
        cf[fineToCoarse[i]] += ff[i];
    }
}


template<typename T>
void matrix::agglomeration::prolongField
(
    Vector<T>       &ff,
    const Vector<T> &cf,
    const label     coarseLevelIndex
) const
{
    const labelField &fineToCoarse = restrictAddressing(coarseLevelIndex);

#ifdef DEBUG
    if (ff.size() != fineToCoarse.size())
    {

        COUT << "Error in prolongField: field does not correspond to level " << coarseLevelIndex
             << " sizes: field = " << ff.size()
             << " level = " << fineToCoarse.size()
             << ENDL;
        ERROR_EXIT;
    }
#endif

    label len = ff.size();

    forAll(i, len)
    {
        ff[i] = cf[fineToCoarse[i]];
    }
}


template<typename T>
void matrix::agglomeration::restrictFaceField
(
    Vector<T>       &cf,
    const Vector<T> &ff,
    const label     fineLevelIndex
) const
{
    const labelField &fineToCoarse = faceRestrictAddressing(fineLevelIndex);

#ifdef DEBUG
    if (ff.size() != fineToCoarse.size())
    {

        COUT << "Error in restrictFaceField: field does not correspond to level " << fineLevelIndex
             << " sizes: field = " << ff.size()
             << " level = " << fineToCoarse.size()
             << ENDL;
        ERROR_EXIT;
    }
#endif

    cf.SET_zero();
    label len = ff.size();

    forAll(i, len)
    {
        label cFace = fineToCoarse[i];

        if (cFace >= 0)
        {
            cf[cFace] += ff[i];
        }
    }
}

} //- namespace UNAP
#endif //- UNAPMATRIX_HPP
