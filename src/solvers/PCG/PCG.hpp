#ifndef PCG_HPP
#define PCG_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class PCG : public matrix::solver
{
private:
  //- use if the CG solver needs to create a null preconditioner
  bool deletePrecondPtr_;

  //- the preconditioner
  matrix::preconditioner *precondPtr_;

public:
  //- constructors
  PCG();

  PCG(matrix::preconditioner &precond);

  //- destructor
  virtual ~PCG()
  {
    if (deletePrecondPtr_)
    {
      delete precondPtr_;
      precondPtr_ = NULL;
      deletePrecondPtr_ = false;
    }
  }

  //- solve the matrix with this solver
  virtual matrix::solverPerformance solve(scalarVector &x,
                                          const matrix &A,
                                          const scalarVector &b) const;
};
}  // namespace UNAP

#endif  //- PCG_HPP
