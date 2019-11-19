#include "printUNAP.hpp"
#include <vector>

using namespace std;

#define output_precision 15

void UNAP::printLDUMatrix(const lduMatrix& A, const char* fileName)
{
	const label nFaces = A.upper().size();
	const label nCells = A.diag().size();
	const scalar* upperPtr = A.upper().begin();
    const scalar* lowerPtr = A.lower().begin();
    const label*  uPtr     = A.upperAddr().begin();
    const label*  lPtr     = A.lowerAddr().begin();
    const scalar* diagPtr  = A.diag().begin();

    label nnz = nCells + nFaces;
    label symm = 0;
    if(A.symm())
    {
    	symm = 1;
    }
    else
    {
    	symm = 0;
    	nnz += nFaces;
    }

    std::ofstream fout;
    FILEOPEN(fout, fileName);
    fout << "nCells: " << nCells << " " << "nnz: " << nnz << " " << "symm: " << symm << std::endl;

    //- diagonal
    forAll(i, nCells)
    {
        fout << i << " " << i << " " << setiosflags(ios::scientific) << setprecision(output_precision) << diagPtr[i] << std::endl;
    }

    //- upper
    forAll(i, nFaces)
    {
        fout << lPtr[i] << " " << uPtr[i] << " "
             << setiosflags(ios::scientific) << setprecision(output_precision)
             << upperPtr[i] << " ";

        //- if there is lower
        if(!symm)
        {
            fout << setiosflags(ios::scientific) << setprecision(output_precision) << lowerPtr[i];
        }

        fout << std::endl;
    }

    FILECLOSE(fout, fileName);

}


void UNAP::printInterfaces(const lduMatrix& A, const char* fileName)
{
	const interfaces& matrixInter = A.matrixInterfaces();
	label numInterfaces = matrixInter.size();

	std::ofstream fout;
	FILEOPEN(fout, fileName);

	fout << "neighbor processors: " << numInterfaces << std::endl;

    for(int i=0; i<numInterfaces; ++i)
    {
        patch& patchI = matrixInter.patchList(i);
        label locSize = patchI.size();
        fout << "processor: " << patchI.neighbProcNo() << " "
             << "size: " << locSize << std::endl;

        for(label facei=0; facei<locSize; facei++)
        {
            fout << patchI.faceCells().begin()[facei] << " "
                 << setiosflags(ios::scientific) << setprecision(output_precision)
                 << patchI.patchCoeffs(facei) << std::endl;
        }
    }

	FILECLOSE(fout, fileName);
}
