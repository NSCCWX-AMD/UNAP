#include "matrixConversion.hpp"

UNAP::lduMatrix& UNAP::coo2ldu
(
	const scalar* dataPtr,
	const label*  rowsPtr,
	const label*  columnPtr,
	const label   nCells,
	const label   size,
	const bool    symm
)
{
	//- input coo matrix must have been arranged in order of row
	//- in each row, data must have been arranged in order of column
	//- structure of matrix must be symmetric
	//- otherwise this part will not work and the results should be wrong

	//- number of no-zero in upper
	const label nZeros = (size - nCells) / 2;

	scalarField diag(nCells);
	scalarField upper(nZeros);

	labelField upperAddr(nZeros);
	labelField lowerAddr(nZeros);

	//- asymmetric case
	scalarField* lowerPtr = NULL;
	labelField*  lowerColPtr = NULL;
	labelField*  upperNbrsPtr = NULL;

	if(!symm)
	{
		lowerPtr = new scalarField(nZeros);
		lowerColPtr = new labelField(nZeros);
		upperNbrsPtr = new labelField(nCells, 0);
	}
	else
	{
		lowerPtr = &upper;
	}

	scalarField& lower = *lowerPtr;
	labelField& lowerCol = *lowerColPtr;
	labelField& upperNbrs = *upperNbrsPtr;


	label upperCount = 0;
	label lowerCount = 0;

	forAll(i, size)
	{
		label row = rowsPtr[i];
		label col = columnPtr[i];

		if(row == col)
		{
			diag[row] = dataPtr[i];
		}
		else if(col > row)
		{
			upperAddr[upperCount] = col;
			lowerAddr[upperCount] = row;
			upper    [upperCount] = dataPtr[i];
			upperCount++;

			if(!symm)
			{
				upperNbrs[row]++;
			}
		}
		else if((col < row) && (!symm))
		{
			lowerCol[lowerCount] = col;
			lower   [lowerCount] = dataPtr[i];
			lowerCount++;
		}
	}

	//- sort lower data
	if(!symm)
	{
		sortData(lower, lowerCol, upperNbrs);
	}

#ifdef DEBUG
	if(upperCount != nZeros)
	{
		UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__
			 << ": the input COO matrix is not a structural symmetric matrix!" << ENDL;
		ERROR_EXIT;
	}
#endif

	lduMatrix* lduAPtr = new lduMatrix
	(
		nCells,
		lowerAddr,
		upperAddr,
		lower,
		diag,
		upper
	);

	if(!symm)
	{
		DELETE_OBJECT_POINTER(lowerPtr)
		DELETE_OBJECT_POINTER(lowerColPtr)
		DELETE_OBJECT_POINTER(upperNbrsPtr)
	}

	return *lduAPtr;
}


UNAP::lduMatrix& UNAP::csr2ldu
(
	const scalar* dataPtr,
	const label*  compRowsPtr,
	const label*  columnPtr,
	const label   nCells,
	const label   size,
	const bool    symm
)
{
	label* rowsPtr = new label[size];

	forAll(i, nCells)
	{
		label nnzInRow = compRowsPtr[i+1] - compRowsPtr[i];
		forAll(j, nnzInRow)
		{
			rowsPtr[compRowsPtr[i] + j] = i;
		}
	}

	lduMatrix& lduA = coo2ldu
	(
		dataPtr,
		rowsPtr,
		columnPtr,
		nCells,
		size,
		symm
	);
	delete [] rowsPtr;
	return lduA;
}


