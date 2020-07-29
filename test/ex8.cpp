#include "matrixConversion.hpp"

//- test for matrix conversion: coo2ldu

//- coo matrix
//- 	1.1    1.2    0.0    1.4    1.5
//-		2.1    2.2    2.3    0.0    0.0
//-     0.0    3.2    3.3    0.0    3.5
//-     4.1    0.0    0.0    4.4    4.5
//-     5.1    0.0    5.3    5.4    5.5

//- coo
//-    data    row    col
//-    1.1     0      0
//-    1.2     0      1
//-    1.4	   0      3
//-    1.5     0      4
//-    2.1     1      0
//-    2.2     1      1
//-    2.3     1      2
//-    3.2     2      1
//-    3.3     2      2
//-    3.5     2      4
//-    4.1     3      0
//-    4.4     3      3
//-    4.5     3      4
//-    5.1     4      0
//-    5.3     4      2
//-    5.4     4      3
//-    5.5     4      4

//- ldu
//- upperAddr:    (1, 3, 4, 2, 4, 4)
//- lowerAddr:    (0, 0, 0, 1, 2, 3)
//- diag:         (1.1, 2.2, 3.3, 4.4, 5.5)
//- upper:        (1.2, 1.4, 1.5, 2.3, 3.5, 4.5)
//- lower:        (2.1, 4.1, 5.1, 3.2, 5.3, 5.4)

using namespace UNAP;

int main()
{
	//- coo
	int row[17] = {0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4, 4};
	int col[17] = {0, 1, 3, 4, 0, 1, 2, 1, 2, 4, 0, 3, 4, 0, 2, 3, 4};
	double data[17] = {1.1, 1.2,      1.4, 1.5,
					   2.1, 2.2, 2.3,
					        3.2, 3.3,      3.5,
					   4.1,           4.4, 4.5,
					   5.1,      5.3, 5.4, 5.5};

	int nCells = 5;
	int size   = 17;
	int nUpperFaces = (size - nCells) / 2;
	bool symm = false;

	lduMatrix& lduA = coo2ldu(data, row, col, nCells, size, symm);


	scalarVector& diag = lduA.diag();
	scalarVector& upper = lduA.upper();
	scalarVector& lower = lduA.lower();
	labelVector& lowerAddr = lduA.lowerAddr();
	labelVector& upperAddr = lduA.upperAddr();

	UNAPCOUT << "upperAddr: " << ENDL;
	forAll(i, nUpperFaces)
	{
		UNAPCOUT << "At i = " << i << ", upperAddr = " << upperAddr[i] << ENDL;
	}

	UNAPCOUT << "lowerAddr: " << ENDL;
	forAll(i, nUpperFaces)
	{
		UNAPCOUT << "At i = " << i << ", lowerAddr = " << lowerAddr[i] << ENDL;
	}

	UNAPCOUT << "diag: " << ENDL;
	forAll(i, nCells)
	{
		UNAPCOUT << "At i = " << i << ", diag = " << diag[i] << ENDL;
	}

	UNAPCOUT << "upper: " << ENDL;
	forAll(i, nUpperFaces)
	{
		UNAPCOUT << "At i = " << i << ", upper = " << upper[i] << ENDL;
	}

	UNAPCOUT << "lower: " << ENDL;
	forAll(i, nUpperFaces)
	{
		UNAPCOUT << "At i = " << i << ", lower = " << lower[i] << ENDL;
	}

	return 0;
}

