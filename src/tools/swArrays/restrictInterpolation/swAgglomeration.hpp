#ifndef SWAgglomeration_HPP
#define SWAgglomeration_HPP

#include "unapMatrix.hpp"
#include "swRestInterStruct.h"

namespace UNAP
{
class swRestInterMap
{
private:
	const matrix::agglomeration& aggl_;
public:
	friend class matrix::agglomeration;

	swRestInterMap(const matrix::agglomeration& aggl)
	:
		aggl_(aggl)
	{}


	static label bandSize_;
	static label minCellsUsingSW_;

	//- restrict
	static bool* restFirstUse_;
	static restStruct* restStructLevels_;

	//- face restrict
	static bool* faceRestFirstUse_;
	static restStruct* faceRestStructLevels_;

	//- interpolate
	static bool* interFirstUse_;
	static interStruct* interStructLevels_;

	//- agglomerate matrix upper
	static bool* aggMatrixUpperFirstUse_;
	static aggMatrixUpperStruct* aggMatrixUpperStructLevels_;

	void initRestInterSize();

	void initRestStruct
	(
		Vector<scalar>& cf,
		const Vector<scalar>& ff,
    	const label fineLevelIndex
	);

	void initInterStruct
	(
		Vector<scalar>& ff,
	    const Vector<scalar>& cf,
	    const label levelIndex
	);

	void agglomerateMatrixUpper
	(
		scalarVector& coarseUpper,
		scalarVector& coarseDiag,
		const scalarVector& fineUpper,
		const label fineLevelIndex
	);

	void initAggMatrixUpperStruct
	(
		scalarVector& coarseUpper,
		scalarVector& coarseDiag,
		const scalarVector& fineUpper,
		const label& fineLevelIndex
	);
};

} //- end namespace UNAP

#endif
