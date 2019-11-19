#ifndef LDUDIAGPRECOND_HPP
#define LDUDIAGPRECOND_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;

class lduDiagPrecond
:
	public matrix::preconditioner
{
private:

	//- the reciprocal diagonal
	scalarField rD_;

public:

	//- constructor
	lduDiagPrecond
	(
		const lduMatrix& A
	);

	//- destructor
	virtual ~lduDiagPrecond()
    {}

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

#endif //- LDUDIAGPRECOND_HPP
