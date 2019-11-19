#ifndef LDUDILUPRECOND_HPP
#define LDUDILUPRECOND_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;

class lduDILUPrecond
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
	lduDILUPrecond
	(
		const lduMatrix& A
	);

	//- destructor
	virtual ~lduDILUPrecond()
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

#endif //- LDUDILUPRECOND_HPP
