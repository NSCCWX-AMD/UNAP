#ifndef LDUDICPRECOND_HPP
#define LDUDICPRECOND_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;

class lduDICPrecond
:
	public matrix::preconditioner
{
private:

	//- the reciprocal diagonal
	scalarField rD_;

	//- lduMatrix
	const lduMatrix* APtr_;

public:

	//- constructor
	lduDICPrecond
	(
		const lduMatrix& A
	);

	//- destructor
	virtual ~lduDICPrecond()
    {}

    //- calculate the reciprocal of the preconditioned diagonal
    static void calcReciprocalD
    (
    	scalarField& rD,
    	const lduMatrix& A
    );

    virtual void precondition
    (
        scalarField& w,
        const scalarField& r
    ) const;

    virtual const scalarField& rD() const
    {
    	return rD_;
    }

};

} //- end namespace UNAP

#endif //- LDUDICPRECOND_HPP
