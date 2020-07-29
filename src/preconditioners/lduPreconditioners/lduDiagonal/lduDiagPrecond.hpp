#ifndef LDUDIAGPRECOND_HPP
#define LDUDIAGPRECOND_HPP

#include "unapMatrix.hpp"

namespace UNAP
{
class lduMatrix;

class lduDiagPrecond : public matrix::preconditioner
{
private:
  //- the reciprocal diagonal
  scalarVector rD_;

public:
  //- constructor
  lduDiagPrecond(const lduMatrix &A);

  //- destructor
  virtual ~lduDiagPrecond() {}

  virtual void precondition(scalarVector &w, const scalarVector &r) const;

  virtual const scalarVector &rD() const { return rD_; }
};

}  // namespace UNAP

#endif  //- LDUDIAGPRECOND_HPP
