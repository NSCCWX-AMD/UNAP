#include "lduMatrix.hpp"

//- test for lduMatrix

using namespace UNAP;

int main()
{
	label nCells = 10;
	scalarField x(nCells, 0.0);

	forAll(i, nCells)
	{
		x[i] = ((scalar)i + 0.2)/0.4;
		UNAPCOUT << "x = " << x[i] << ENDL;
	}

	labelField lowerAddr(nCells);
	labelField upperAddr(nCells);

	forAll(i, nCells)
	{
		lowerAddr[i] = i == nCells-1? 0       : i;
		upperAddr[i] = i == nCells-1? nCells-1 : i+1;
	}

	scalarField lower(nCells, -1.0);
	scalarField upper(nCells, -1.0);
	scalarField diag (nCells, 4.0);

	lduMatrix lduM(nCells,
		   		   lowerAddr,
		   		   upperAddr,
		   		   lower,
		   		   diag,
		   		   upper);

	matrix *A = &lduM;

	scalarField Apsi(nCells, 0.0);
    A->spMV(Apsi, x);

    forAll(i, nCells)
    {
    	UNAPCOUT << "i = " << i << ", Apsi = " << Apsi[i] << ENDL;
    }

    //- for checking the correct
}