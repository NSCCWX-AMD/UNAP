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
		scalarVector       &x,
        const matrix      &A,
        const scalarVector &b,
        const label       nSweeps
	) const;

	void smooth
	(
		scalarVector       &x,
        const lduMatrix   &A,
        const scalarVector &b,
        const label       nSweeps
	) const;

	virtual void init() const
	{}

};

} //- end namespace UNAP

#endif //- LDUGAUSSSEIDELSMOOTHER_HPP
