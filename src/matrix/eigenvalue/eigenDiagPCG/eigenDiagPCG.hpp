#ifndef EIGENDIAGPCG_HPP
#define EIGENDIAGPCG_HPP

#include "unapMatrix.hpp"
#include <vector>

namespace UNAP
{
class eigenDiagPCG
{
private:
	//- the maximum eigenvalue
	mutable scalar maxEigenValue_;

	//- Sturm Sequence Method
    //- find polynomial vector p(i) to find eigen value
    void h14Sturm
    (
        scalar** TriMatrix,  	/// tridiagonal matrix
        const scalar lamb,   	/// guess eigen value (used in bisection method)
        std::vector<scalar>& p, /// a pvector polynomials (pi_1 to pi_n)
        label& s,            	/// the number of sign change
        const label nPCGs    	/// nDiagPCGs_
    ) const;

    //- form the matrix from alphas and betas, used to calculate eigenvalue
    //- algorithm see Ref.[3]
    void computeValueForMatrix
    (
        const scalarVector& alphas,  /// alphas from nPCGs times of PCG
        const scalarVector& betas,   /// betas from nPCGs times of PCG
        scalar** TriMatrix,         /// tridiagonal matrix
        const label nPCGs           /// nDiagPCGs_
    ) const;

    //- estimate the range of eigenvalues
    void determineEigRange
    (
        scalar** TriMatrix,   /// tridiagonal matrix
        scalar& xBegin,       /// left bound of eigenvalues
        scalar& xEnd,         /// right bound of eigenvalues
        const label nPCGs     /// nDiagPCGs_
    ) const;

    //- compute eigenvalue
    void computeMaxEig
    (
        const scalarVector& alphas,   /// alphas from nPCGs times of PCG
        const scalarVector& betas,    /// betas from nPCGs times of PCG
        const label nPCGs,           /// nDiagPCGs_
        const label k               /// kth largest eigenvalue, default k=1 for the maximum eigenvalue
    ) const;

    //- diagPCG loops
    void diagPCGLoops
    (
    	const matrix& A,
		scalarVector& x,
		const scalarVector& b,
		const matrix::preconditioner& precond,
		const label nDiagPCGs,
		scalarVector& alphas,
		scalarVector& betas
    ) const;

    template <typename T> T** allocateSym2D(label nPCGs) const;
    template <typename T> void deleteSym2D(T** arr, label nPCGs) const;

public:
	eigenDiagPCG
	(
		const matrix& A,
		scalarVector& x,
		const scalarVector& b,
		const matrix::preconditioner& precond,
		const label nDiagPCGs
	);

	//- return the maximum eigenvalue
	scalar maxEigenValue() const
	{
		return maxEigenValue_;
	}
};

// allocate memory
template <typename T>
T** UNAP::eigenDiagPCG::allocateSym2D (label n) const
{
    T** arr2D;
    arr2D = new T* [n];

    for(label i=0; i<n; i++)
    {
        arr2D[i] = new T[n];
    }

    return (T**)arr2D;
}


// free memory
template <typename T>
void UNAP::eigenDiagPCG::deleteSym2D(T** arr, label n) const
{
    if(arr != NULL)
    {
        for(label i=0; i<n; i++)
        {
            delete [] arr[i];
        }

        delete [] arr;
        arr = NULL;
    }
}


} //- namespace UNAP
#endif //-EIGENDIAGCG_HPP
