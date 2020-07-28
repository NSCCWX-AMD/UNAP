#include "lduAgglomeration.hpp"

void UNAP::lduAgglomeration::agglomerateMatrix
(
	const label fineLevelIndex
)
{
	//- get fine matrix
	const lduMatrix& fineMatrix = matrixLevel(fineLevelIndex);

	//- set the coarse level matrix
	//- here should give exact lduMatrix parameters(upperAddr, lowerAddr, ...),
	//  except for coefficients
	lduMatrix& coarseMatrix = coarseMatrixLevels_[fineLevelIndex];

	//- get face restriction map for current level
	const labelField& faceRestrictAddr = faceRestrictAddressing(fineLevelIndex);

	//- coarse matrix diagonal initialized by restricting the finer mesh diagonal
    coarseMatrix.SET_diag(nCells_[fineLevelIndex]);

	scalarField& coarseDiag = coarseMatrix.diag();

    restrictField(coarseDiag, fineMatrix.diag(), fineLevelIndex);

    //- check if matrix is asymmetric and if so agglomerate both upper and lower
    //  coefficients ...
    if(!fineMatrix.symm())
    {
        //- get off-diagonal matrix coefficients
        const scalarField& fineUpper = fineMatrix.upper();
        const scalarField& fineLower = fineMatrix.lower();

        //- coarse matrix upper coefficients
        coarseMatrix.SET_upper(coarseMatrix.upperAddr().size());
        coarseMatrix.SET_lower(coarseMatrix.lowerAddr().size());
        scalarField& coarseUpper = coarseMatrix.upper();
        scalarField& coarseLower = coarseMatrix.lower();

        const labelField& restrictAddr = restrictAddressing(fineLevelIndex);

        const labelField& l  = fineMatrix.lowerAddr();
        const labelField& cl = coarseMatrix.lowerAddr();
        const labelField& cu = coarseMatrix.upperAddr();

        forAll(fineFacei, faceRestrictAddr.size())
        {
            label cFace = faceRestrictAddr[fineFacei];

            if (cFace >= 0)
            {
                //- check the orientation of the fine-face relative to the
                //  coarse face it is being agglomerated into
                if (cl[cFace] == restrictAddr[l[fineFacei]])
                {
                    coarseUpper[cFace] += fineUpper[fineFacei];
                    coarseLower[cFace] += fineLower[fineFacei];
                }
                else if (cu[cFace] == restrictAddr[l[fineFacei]])
                {
                    coarseUpper[cFace] += fineLower[fineFacei];
                    coarseLower[cFace] += fineUpper[fineFacei];
                }
                else
                {
                    UNAPCOUT << "Error in agglomerateMatrix: Inconsistent addressing between "
                    	 << "fine and coarse grids." << ENDL;
                    ERROR_EXIT;
                }
            }
            else
            {
                //- add the fine face coefficients into the diagonal.
                coarseDiag[-1 - cFace] +=
                    fineUpper[fineFacei] + fineLower[fineFacei];
            }
        }
    }
    else //- ... Otherwise it is symmetric so agglomerate just the upper
    {
        //- get off-diagonal matrix coefficients
        const scalarField& fineUpper = fineMatrix.upper();

        //- coarse matrix upper coefficients
        coarseMatrix.SET_upper(coarseMatrix.upperAddr().size());
        scalarField& coarseUpper = coarseMatrix.upper();

        coarseMatrix.setSymm();

        forAll(fineFacei, faceRestrictAddr.size())
        {
            label cFace = faceRestrictAddr[fineFacei];

            if (cFace >= 0)
            {
                coarseUpper[cFace] += fineUpper[fineFacei];
            }
            else
            {
                //- add the fine face coefficient into the diagonal.
                coarseDiag[-1 - cFace] += 2*fineUpper[fineFacei];
            }
        }
    }

    if(PARRUN)
    {
        const interfaces& fineInterfaces = fineMatrix.matrixInterfaces();
        interfaces& coarseInterfaces = coarseMatrix.matrixInterfaces();

        const label interfacesSize = fineInterfaces.size();
        forAll(inti, interfacesSize)
        {
            const scalarField& finePatchCoeffs = fineInterfaces.patchList(inti).patchCoeffs();
            const label coarsePatchFacesSize = coarseInterfaces.patchList(inti).size();
            const label finePatchFacesSize = fineInterfaces.patchList(inti).size();

            scalarField* coarsePatchCoeffsPtr = new scalarField(coarsePatchFacesSize);
            scalarField& coarsePatchCoeffs = *coarsePatchCoeffsPtr;

            labelField& faceRestrictAddressing = coarseInterfaces.patchList(inti).faceRestrictAddressing();

            forAll(ffi, finePatchFacesSize)
            {
                coarsePatchCoeffs[faceRestrictAddressing[ffi]] += finePatchCoeffs[ffi];
            }

            coarseInterfaces.patchList(inti).patchCoeffs(coarsePatchCoeffs);
        }
    }
}
