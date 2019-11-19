#ifndef LDUGAUSSSEIDELSMOOTHER_HPP
#define LDUGAUSSSEIDELSMOOTHER_HPP

#include "lduMatrix.hpp"

namespace UNAP
{

class lduGaussSeidelSmoother
:
	public matrix::smoother
{

private:


public:

	//- constructors
	lduGaussSeidelSmoother()
	{
		// if(!MYID)
  //       {
  //           COUT << "Gauss-Seidel smoother used!" << ENDL;
  //       }
	}

	//- destructor
	virtual ~lduGaussSeidelSmoother()
	{}

	//- smooth the solution for a given number of sweeps
	virtual void smooth
	(
		scalarField       &x,
        const matrix      &A,
        const scalarField &b,
        const label       nSweeps
	) const;

	void smooth
	(
		scalarField       &x,
        const lduMatrix   &A,
        const scalarField &b,
        const label       nSweeps
	) const;

	virtual void init() const
	{}

};

} //- end namespace UNAP

#endif //- LDUGAUSSSEIDELSMOOTHER_HPP
