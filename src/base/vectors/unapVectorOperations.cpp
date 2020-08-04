#include "unapVector.hpp"

namespace UNAP
{
scalar dot(const scalarVector &v1, const scalarVector &v2)
{
#ifdef DEBUG
  if (v1.size() != v2.size())
  {
    v1.getCommunicator()->log()
        << "ERROR in " << __FILE__ << " " << __LINE__
        << ": The length of two vectors is not same!" << ENDL;
    ERROR_EXIT;
  }

  if (v1.getCommunicator() != v2.getCommunicator())
    v1.getCommunicator()->log
        << "ERROR in " << __FILE__ << " " << __LINE__
        << ": The Communicator of two vectors is not same!" << ENDL;

#endif

  scalar res = 0.0;
  const label len = v1.size();
  const scalar *val1 = v1.values();
  const scalar *val2 = v2.values();

  forAll(i, len) { res += val1[i] * val2[i]; }

  reduceSum(&res, v1.getCommunicator());
  return res;
}

}  // namespace UNAP
