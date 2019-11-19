#include "matrixConversion.hpp"
#include "printUNAP.hpp"

void UNAP::sortData
(
	scalarField& data,
	const labelField&  order,
	const labelField&  cellFaces
)
{
	const label nCells = cellFaces.size();

	labelField cellFaceOffsets(nCells + 1);

	cellFaceOffsets[0] = 0;
	forAll(celli, nCells)
	{
		cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + cellFaces[celli];
	}

	//- size of faces in upper
	const label nFaces = data.size();

	scalarField dataOld(data);

	//- count the data input
	labelField count(nCells, 0);

	forAll(facei, nFaces)
	{
		label cellLoc = order[facei];
		data[cellFaceOffsets[cellLoc]+count[cellLoc]] = dataOld[facei];
		count[cellLoc]++;
	}
}


void UNAP::reorderCOO
(
	scalar* dataPtr,
	label*  rowsPtr,
	label*  colsPtr,
	const label nCells,
	const label size
)
{
	labelField countsInRow(nCells);

	forAll(i, size)
	{
		countsInRow[rowsPtr[i]]++;
	}

	labelField countsInRowOffsets(nCells+1);

	forAll(i, nCells)
	{
		countsInRowOffsets[i+1] = countsInRowOffsets[i] + countsInRow[i];
	}

	labelField posInRow(size);
	forAll(i, size)
	{
		label rowLocal = rowsPtr[i];
		label rowStartPos = countsInRowOffsets[rowLocal];
		label localRowSize = countsInRowOffsets[rowLocal+1] - countsInRowOffsets[rowLocal];

		forAll(k, localRowSize)
		{
			if(colsPtr[i] > colsPtr[rowStartPos+k])
			{
				posInRow[i]++;
			}
		}
	}

	labelField rowTemp(size);
	labelField colTemp(size);
	scalarField valTemp(size);

	forAll(i, size)
	{
		label rowLocal = rowsPtr[i];
		label rowStartPos = countsInRowOffsets[rowLocal];
		label pos = rowStartPos + posInRow[i];

		valTemp[pos] = dataPtr[i];
		rowTemp[pos] = rowsPtr[i];
		colTemp[pos] = colsPtr[i];
	}

	forAll(i, size)
	{
		dataPtr[i] = valTemp[i];
		rowsPtr[i] = rowTemp[i];
		colsPtr[i] = colTemp[i];
	}
}


void UNAP::reorderValue
(
	scalar* val,
	const label* newOrder,
	const label  size
)
{
	scalarField valTemp(val, size, false);

	forAll(i, size)
	{
		val[newOrder[i]] = valTemp[i];
	}
}


void UNAP::reorderUFace
(
	label*  rowsPtr,
	label*  colsPtr,
	const label   nCells,
	const label   size,
	label*  newOrder
)
{
	labelField countsInRow(nCells);

	forAll(i, size)
	{
		countsInRow[rowsPtr[i]]++;
	}

	labelField countsInRowOffsets(nCells+1);

	forAll(i, nCells)
	{
		countsInRowOffsets[i+1] = countsInRowOffsets[i] + countsInRow[i];
	}

	labelField posInRow(size);
	forAll(i, size)
	{
		label rowLocal = rowsPtr[i];
		label rowStartPos = countsInRowOffsets[rowLocal];
		label localRowSize = countsInRowOffsets[rowLocal+1] - countsInRowOffsets[rowLocal];

		forAll(k, localRowSize)
		{
			if(colsPtr[i] > colsPtr[rowStartPos+k])
			{
				posInRow[i]++;
			}
		}
	}

	labelField rowTemp(size);
	labelField colTemp(size);

	forAll(i, size)
	{
		label rowLocal = rowsPtr[i];
		label rowStartPos = countsInRowOffsets[rowLocal];
		label pos = rowStartPos + posInRow[i];

		newOrder[i] = pos;
		rowTemp[pos] = rowsPtr[i];
		colTemp[pos] = colsPtr[i];
	}

	forAll(i, size)
	{
		rowsPtr[i] = rowTemp[i];
		colsPtr[i] = colTemp[i];
	}
}


void UNAP::reorderLFace
(
	label*  rowsPtr,
	label*  colsPtr,
	const label   nCells,
	const label   size,
	label*  newOrder
)
{
	labelField countsInCol(nCells);

	forAll(i, size)
	{
		countsInCol[colsPtr[i]]++;
	}

	labelField countsInColOffsets(nCells+1);

	forAll(i, nCells)
	{
		countsInColOffsets[i+1] = countsInColOffsets[i] + countsInCol[i];
		countsInCol[i] = 0;
	}

	labelField rowTemp(size);
	labelField colTemp(size);

	forAll(i, size)
	{
		label colLocal = colsPtr[i];
		label colStartPos = countsInColOffsets[colLocal];
		label pos = colStartPos + countsInCol[colLocal];

		newOrder[i]  = pos;
		rowTemp[pos] = rowsPtr[i];
		colTemp[pos] = colsPtr[i];
		countsInCol[colLocal]++;
	}

	forAll(i, size)
	{
		rowsPtr[i] = rowTemp[i];
		colsPtr[i] = colTemp[i];
	}
}

