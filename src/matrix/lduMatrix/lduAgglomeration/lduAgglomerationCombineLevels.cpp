#include "lduAgglomeration.hpp"

void UNAP::lduAgglomeration::combineLevels(const label curLevel)
{
	label prevLevel = curLevel - 1;

	//- set the previous level nCells to the current
    nCells_[prevLevel] = nCells_[curLevel];

    //- map the restrictAddressing from the coarser level into the previous
    //  finer level

    const labelField& curResAddr = restrictAddressing_[curLevel];
    labelField& prevResAddr = restrictAddressing_[prevLevel];

    const labelField& curFaceResAddr = faceRestrictAddressing_[curLevel];
    labelField& prevFaceResAddr = faceRestrictAddressing_[prevLevel];

    forAll(i, prevFaceResAddr.size())
    {
        if (prevFaceResAddr[i] >= 0)
        {
            prevFaceResAddr[i] = curFaceResAddr[prevFaceResAddr[i]];
        }
        else
        {
            prevFaceResAddr[i] = -curResAddr[-prevFaceResAddr[i] - 1] - 1;
        }
    }

    //- delete the restrictAddressing for the coarser level
    faceRestrictAddressing_.removeLevel(curLevel);

    forAll(i, prevResAddr.size())
    {
        prevResAddr[i] = curResAddr[prevResAddr[i]];
    }

    //- delete the restrictAddressing for the coarser level
    restrictAddressing_.removeLevel(curLevel);

    //- delete the matrix addressing and coefficients from the previous level
    //  and replace with the corresponding entries from the coarser level
    coarseMatrixLevels_.removeLevel(prevLevel);



    //- same for the lduInterfaceFields taking care to delete the sub-entries
    //  held on List<T*>
    // const lduInterfacePtrsList& curInterLevel = interfaceLevels_[curLevel+1];
    // lduInterfacePtrsList& prevInterLevel = interfaceLevels_[prevLevel+1];

    // forAll(prevInterLevel, inti)
    // {
    //     if (prevInterLevel.set(inti))
    //     {
    //         refCast<GAMGInterface>(const_cast<lduInterface&>
    //         (
    //             prevInterLevel[inti]
    //         )).combine(refCast<const GAMGInterface>(curInterLevel[inti]));

    //         delete curInterLevel(inti);
    //     }
    // }

    // interfaceLevels_.set(curLevel+1, NULL);
}

