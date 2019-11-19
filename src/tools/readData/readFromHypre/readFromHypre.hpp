#ifndef READFROMHYPRE_HPP
#define READFROMHYPRE_HPP

#include "lduMatrix.hpp"
#include <map>

namespace UNAP
{

void constructLDUMatrixFromHypre(lduMatrix& lduA, const char* fileName);

void sortInterFaces
(
	scalarField& val,
	labelField&  row,
	labelField&  col,
	const labelField& faceToProcNO,
	const label  faceSize,
	const labelField& faceStart,
	std::map<int, int>& mapProcNO,
	const labelField& neiProcNo,
	const label  procSize,
	const labelField& globalRowStart,
	const labelField& globalRowEnd
);

void constructLDUInterfacesFromHypre
(
	lduMatrix&  lduA,
	const label nNeiProcs,
	const labelField&  destRank,
	const labelField&  locPosition,
	const labelField&  faceCells,
	const scalarField& data
);

void constructVectorFromHypre(scalarField& b, const char* fileName);
}
#endif
