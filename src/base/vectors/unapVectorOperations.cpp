#include "unapVector.hpp"

namespace UNAP
{

scalar dot
(
	const scalarVector &v1,
	const scalarVector &v2
)
{
#ifdef DEBUG
	if(v1.size() != v2.size())
	{
		UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__
			 << ": The length of two vectors is not same!" << ENDL;
		ERROR_EXIT;
	}
#endif

	scalar res = 0.0;
	const label len = v1.size();
	const scalar *val1 = v1.values();
	const scalar *val2 = v2.values();


	forAll(i, len)
	{
		res += val1[i] * val2[i];
	}

	reduceSum(res);
	return res;
}

}

