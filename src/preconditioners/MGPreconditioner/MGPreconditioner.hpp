#ifndef MGPRECONDITIONER_HPP
#define MGPRECONDITIONER_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;
class MGSolver;

class MGPrecond : matrix::preconditioner
{
private:
  //- number of V-cycles to perform
  label nVcycles_;

  //- matrix pointer
  const lduMarix *APtr_;

  //- MG solver pointer
  const MGSolver *MGSolverPtr_;

public:
  //- constructor
  MGPrecond(const lduMatrix &A);

  //- destructor
  virtual ~MGPrecond() {}

  //- return wA the preconditioned form of residual rA
  virtual void precondition(scalarVector &wA, const scalarVector &rA) const;

  void set_nVcycles(const label n) { nVcycles_ = n; }
};

}  // namespace UNAP

#endif  //- MGPRECONDITIONER_HPP
