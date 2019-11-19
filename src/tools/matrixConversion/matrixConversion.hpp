#ifndef MATRIXCONVERSION_HPP
#define MATRIXCONVERSION_HPP

#include "lduMatrix.hpp"

namespace UNAP
{

void sortData
(
	scalarField& data,
	const labelField& order,
	const labelField& cellFaces
);


//- reorder COO matrix, to make it ordered in columns in every row
//- to use this, the input COO matrix must be ordered in rows
//- row, col must be started from 0
void reorderCOO
(
	scalar* dataPtr,
	label*  rowsPtr,
	label*  columnPtr,
	const label   nCells,
	const label   size
);


void reorderValue
(
	scalar* val,
	const label* newOrder,
	const label  size
);


void reorderUFace
(
	label*  rowsPtr,
	label*  columnPtr,
	const label   nCells,
	const label   size,
	label*  newOrder
);


void reorderLFace
(
	label*  rowsPtr,
	label*  columnPtr,
	const label   nCells,
	const label   size,
	label*  newOrder
);


//- only for diagonal part
lduMatrix& coo2ldu
(
	const scalar* dataPtr,
	const label*  rowsPtr,
	const label*  columnPtr,
	const label   nCells,
	const label   size,
	const bool    symm  //- symm refers to the data
);

//- only for diagonal part
lduMatrix& csr2ldu
(
	const scalar* dataPtr,
	const label*  compRowsPtr,
	const label*  columnPtr,
	const label   nCells,
	const label   size,
	const bool    symm  //- symm refers to the data
);

} //- end namespace UNAP

#endif //- MATRIXCONVERSION_HPP
