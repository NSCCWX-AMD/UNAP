#ifndef LDUGAUSSSEIDELSMOOTHER_HPP
#define LDUGAUSSSEIDELSMOOTHER_HPP

#include "lduMatrix.hpp"

namespace UNAP
{
class lduGaussSeidelSmoother : public matrix::smoother
{
 private:
 public:
  //- constructors
  lduGaussSeidelSmoother(Communicator *other_comm)
      : matrix::smoother(other_comm)
  {
  }

  //- destructor
  virtual ~lduGaussSeidelSmoother() {}

  //- smooth the solution for a given number of sweeps
  virtual void smooth(scalarVector &x,
                      const matrix &A,
                      const scalarVector &b,
                      const label nSweeps) const;

  void smooth(scalarVector &x,
              const lduMatrix &A,
              const scalarVector &b,
              const label nSweeps) const;

  virtual void init() const {}
};

}  // namespace UNAP

#endif  //- LDUGAUSSSEIDELSMOOTHER_HPP
