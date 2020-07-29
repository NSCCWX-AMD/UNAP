#ifndef CHEBYSMOOTHER_CPP
#define CHEBYSMOOTHER_CPP

#include "lduMatrix.hpp"

namespace UNAP
{
class chebySmoother : public matrix::smoother
{
private:
  //- Number of PCGs to get alphas and betas, used to calculate eigenvalue
  label nDiagPCGs_;

  //- Upper bound of the bounding ellipse of the eigenvalues of the matrix A
  mutable scalar maxEigPCG_;

  //- Detect if eigenvalue has been calculated,
  //- if not, Chebyshev smoother will call PCGs
  mutable bool eigFirstTimeComputed_;

  //- Estimate the minEig by maxEig divided by eigRatio
  //- minEig = maxEig / eigRatioCheby_
  scalar eigRatioCheby_;

  //- factor to enlarge the maxEig, because
  //- the maxEig obtained is usually underestimated
  //- maxEig *= boostFactorCheby_
  scalar boostFactorCheby_;

  //- not used here
  scalar eigRatioCoarest_;

  //-detect if using preSmooth_
  //-if using, first Amul operation can be ignored
  label preSmoothUsing_;

public:
  chebySmoother()
      : nDiagPCGs_(10),
        maxEigPCG_(0.0),
        eigFirstTimeComputed_(true),
        eigRatioCheby_(1.5),     // 30
        boostFactorCheby_(1.1),  // 1.05
        eigRatioCoarest_(1.0),
        preSmoothUsing_(0)
  {
    // if(!MYID)
    // {
    //     COUT << "Chebyshev smoother used!" << ENDL;
    // }
  }

  //- smooth the solution for a given number of sweeps
  virtual void smooth(scalarVector &x,
                      const matrix &A,
                      const scalarVector &b,
                      const label nSweeps) const;

  void smooth(scalarVector &x,
              const lduMatrix &A,
              const scalarVector &b,
              const label nSweeps) const;

  virtual void init() const { eigFirstTimeComputed_ = true; }
};
}  // namespace UNAP

#endif  //- CHEBYSMOOTHER_CPP
