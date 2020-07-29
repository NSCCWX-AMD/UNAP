#ifndef PBICGSTAB_HPP
#define PBICGSTAB_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class PBiCGStab
:
	public matrix::solver
{
private:

	//- use if the CG solver needs to create a null preconditioner
	bool deletePrecondPtr_;

	//- the preconditioner
	matrix::preconditioner* precondPtr_;

public:

	//- constructors
	PBiCGStab();

	PBiCGStab
	(
		matrix::preconditioner& precond
	);

	//- destructor
	virtual ~PBiCGStab()
	{
		if(deletePrecondPtr_)
		{
			delete precondPtr_;
			precondPtr_ = NULL;
			deletePrecondPtr_ = false;
		}
	}

	//- solve the matrix with this solver
	virtual matrix::solverPerformance solve
	(
		scalarVector& x,
		const matrix& A,
		const scalarVector& b
	) const;

	void SET_preconditioner(matrix::preconditioner& precond)
	{
		DELETE_OBJECT_POINTER(precondPtr_);
		precondPtr_ = &precond;
	}

};
} //- namespace UNAP



#endif //- PBICGSTAB_HPP
