#include "lduAgglomeration.hpp"
#include "printUNAP.hpp"

void UNAP::lduAgglomeration::agglomerate
(
    const scalarField& faceWeights
)
{
	//- start geometric agglomeration from the given faceWeights
	scalarField* faceWeightsPtr = const_cast<scalarField*>(&faceWeights);

	//- agglomerate until the required number of cells in the coarsest level
    //  is reached

	label nPairLevels = 0;
    label nCreatedLevels = 0;

    while (nCreatedLevels < maxLevels_ - 1)
    {
    	label nCoarseCells = -1;

        labelField* finalAgglomPtr = &agglomerate
        (
            nCoarseCells,
            matrixLevel(nCreatedLevels),
            *faceWeightsPtr
        );

    	if (continueAgglomerating(nCoarseCells))
        {
            nCells_[nCreatedLevels] = nCoarseCells;
            restrictAddressing_.setLevel(nCreatedLevels, *finalAgglomPtr);
        }
        else
        {
            DELETE_OBJECT_POINTER(finalAgglomPtr)
            break;
        }

    	agglomerateLduAddressing(nCreatedLevels);

    	//- agglomerate the faceWeights field for the next level
    	{
    		label nCoarseFaces = coarseMatrixLevels_[nCreatedLevels].upperAddr().size();

    		scalarField* aggFaceWeightsPtr = new scalarField(nCoarseFaces, 0.0);

    		restrictFaceField
            (
                *aggFaceWeightsPtr,
                *faceWeightsPtr,
                nCreatedLevels
            );

            if (nCreatedLevels)
            {
                delete faceWeightsPtr;
                faceWeightsPtr = NULL;
            }

            faceWeightsPtr = aggFaceWeightsPtr;
    	}

    	if (nPairLevels % mergeLevels_)
        {
            combineLevels(nCreatedLevels);
        }
        else
        {
            nCreatedLevels++;
        }

        nPairLevels++;
    }

    // shrink the storage of the levels to those created
    compactLevels(nCreatedLevels);

    // delete temporary geometry storage
    if (nCreatedLevels)
    {
        delete faceWeightsPtr;
    }

#ifdef DEBUG
    
        UNAPCOUT << nCreatedLevels << " coarse levels created!" << ENDL;
    
#endif
}


UNAP::labelField& UNAP::lduAgglomeration::agglomerate
(
	label&             nCoarseCells,
    const lduMatrix&   fineA,
    const scalarField& faceWeights
)
{
	const label nFineCells = fineA.size();

	const labelField& upperAddr = fineA.upperAddr();
	const labelField& lowerAddr = fineA.lowerAddr();

	//- for each cell calculate faces
	labelField cellFaces(upperAddr.size() + lowerAddr.size());
	labelField cellFaceOffsets(nFineCells + 1);

	//- memory management
	{
		labelField nNbrs(nFineCells, 0);

		//- get number of faces in each row
        forAll(facei, upperAddr.size())
		{
			nNbrs[upperAddr[facei]]++;
		}

		forAll(facei, lowerAddr.size())
		{
			nNbrs[lowerAddr[facei]]++;
		}

		cellFaceOffsets[0] = 0;

		forAll(celli, nNbrs.size())
        {
            cellFaceOffsets[celli+1] = cellFaceOffsets[celli] + nNbrs[celli];
        }

        //- reset the whole list to use as counter
        nNbrs = 0;

        //- faces are stored in cell order, including both upper and lower
        forAll(facei, upperAddr.size())
        {
        	cellFaces
            [
                cellFaceOffsets[upperAddr[facei]] + nNbrs[upperAddr[facei]]
            ] = facei;

            nNbrs[upperAddr[facei]]++;
        }

        forAll(facei, lowerAddr.size())
        {
            cellFaces
            [
                cellFaceOffsets[lowerAddr[facei]] + nNbrs[lowerAddr[facei]]
            ] = facei;

            nNbrs[lowerAddr[facei]]++;
        }
	}

	//- go through the faces and create clusters
    labelField* coarseCellMapPtr = new labelField(nFineCells, -1);
    labelField& coarseCellMap =  *coarseCellMapPtr;

	nCoarseCells = 0;

	for(label celli=0; celli<nFineCells; celli++)
	{
		//- change cell ordering depending on direction for this level
        celli = forward_ ? celli : nFineCells - celli - 1;

        if (coarseCellMap[celli] < 0)
		{
			label matchFaceNo = -1;
			scalar maxFaceWeight = SMALL;

			//- check faces to find ungrouped neighbor with largest face weight
			for
            (
                label faceOs = cellFaceOffsets[celli];
                faceOs < cellFaceOffsets[celli+1];
                faceOs++
            )
            {
                label facei = cellFaces[faceOs];

                //- I don't know whether the current cell is owner or neighbor.
                //- Therefore I'll check both sides
                if
                (
                    coarseCellMap[upperAddr[facei]] < 0
                 && coarseCellMap[lowerAddr[facei]] < 0
                 && faceWeights[facei] > maxFaceWeight
                )
                {
                    //- match found. Pick up all the necessary data
                    matchFaceNo = facei;
                    maxFaceWeight = faceWeights[facei];
                }
            }

            if (matchFaceNo >= 0)
            {
                //- make a new group
                coarseCellMap[upperAddr[matchFaceNo]] = nCoarseCells;
                coarseCellMap[lowerAddr[matchFaceNo]] = nCoarseCells;
                nCoarseCells++;
            }
            else
            {
            	//- no match. Find the best neighboring cluster and
                //  put the cell there
                label clusterMatchFaceNo = -1;
                scalar clusterMaxFaceCoeff = SMALL;

                for
                (
                    label faceOs=cellFaceOffsets[celli];
                    faceOs<cellFaceOffsets[celli+1];
                    faceOs++
                )
                {
                    label facei = cellFaces[faceOs];

                    if (faceWeights[facei] > clusterMaxFaceCoeff)
                    {
                        clusterMatchFaceNo = facei;
                        clusterMaxFaceCoeff = faceWeights[facei];
                    }
                }

                if (clusterMatchFaceNo >= 0)
                {
                    //- add the cell to the best cluster
                    coarseCellMap[celli] = MAX
                    (
                        coarseCellMap[upperAddr[clusterMatchFaceNo]],
                        coarseCellMap[lowerAddr[clusterMatchFaceNo]]
                    );
                }
            }
		}
	}


	//- check that all cells are part of clusters,
    //- if not create single-cell "clusters" for each
    for (label celli=0; celli<nFineCells; celli++)
    {
        //- change cell ordering depending on direction for this level
        celli = forward_ ? celli : nFineCells - celli - 1;

        if (coarseCellMap[celli] < 0)
        {
            coarseCellMap[celli] = nCoarseCells;
            nCoarseCells++;
        }
    }

    //- reverse the map ordering to improve the next level of agglomeration
    //  (doesn't always help and is sometimes detrimental)
    if (!forward_)
    {
        nCoarseCells--;

        forAll(celli, coarseCellMap.size())
        {
            coarseCellMap[celli] = nCoarseCells - coarseCellMap[celli];
        }

        nCoarseCells++;
    }

    return coarseCellMap;
}
