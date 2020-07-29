#ifndef LDUDICPRECOND_HPP
#define LDUDICPRECOND_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;

class lduDICPrecond : public matrix::preconditioner
{
private:
  //- the reciprocal diagonal
  scalarVector rD_;

  //- lduMatrix
  const lduMatrix *APtr_;

public:
  //- constructor
  lduDICPrecond(const lduMatrix &A);

  //- destructor
  virtual ~lduDICPrecond() {}

  //- calculate the reciprocal of the preconditioned diagonal
  static void calcReciprocalD(scalarVector &rD, const lduMatrix &A);

  virtual void precondition(scalarVector &w, const scalarVector &r) const;

  virtual const scalarVector &rD() const { return rD_; }
};

}  // namespace UNAP

#endif  //- LDUDICPRECOND_HPP
