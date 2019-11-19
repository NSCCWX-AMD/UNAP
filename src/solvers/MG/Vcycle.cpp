#include "MG.hpp"

#ifdef SW_SLAVE
#include "vectorOps.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

void UNAP::MGSolver::initVcycle
(
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources
) const
{
    forAll(leveli, agglomeration_.size())
    {
        label nCells = agglomeration_.coarseMatrixLevels(leveli).size();
        scalarField* coarseCorrFieldsLeveliPtr = new scalarField(nCells);
        scalarField* coarseSourcesLeveliPtr = new scalarField(nCells);
        coarseCorrFields.setLevel(leveli, *coarseCorrFieldsLeveliPtr);
        coarseSources.setLevel(leveli, *coarseSourcesLeveliPtr);
    }
}

void UNAP::MGSolver::Vcycle
(
    scalarField& psi,
    const scalarField& source,
    scalarField& Apsi,
    scalarField& finestCorrection,
    scalarField& finestResidual,
    PtrList<scalarField>& coarseCorrFields,
    PtrList<scalarField>& coarseSources
) const
{
#ifdef SW_SLAVE
    MVM_Arrays arrays1;
#endif
    const label coarsestLevel = agglomeration_.size() - 1;

    //- restrict finest grid residual for the next level up
    if(coarsestLevel >= 0)
        agglomeration_.restrictField(coarseSources[0], finestResidual, 0);

    //- residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        //- if the optional pre-smoothing sweeps are selected
        //  smooth the coarse-grid field for the restricted source

        label nCells = coarseCorrFields[leveli].size();
        scalar* coarseCorrFieldsPtr = coarseCorrFields[leveli].values();

        if (nPreSweeps_)
        {
            coarseCorrFields[leveli] = 0.0;

            smoothers_[leveli+1].smooth
            (
                coarseCorrFields[leveli],
                agglomeration_.coarseMatrixLevels(leveli),
                coarseSources[leveli],
                nPreSweeps_ + leveli
            );

            scalarField ACf(nCells);

            //- scale coarse-grid correction field
            //  but not on the coarsest level because it evaluates to 1
            if (scaleCorrection_ && leveli < coarsestLevel - 1)
            {
                scalar sf = scalingFactor
                (
                    ACf,
                    coarseCorrFields[leveli],
                    coarseSources[leveli],
                    agglomeration_.coarseMatrixLevels(leveli)
                );

                IFNOT_SWACC
                {
                    forAll(cellI, nCells)
                    {
                        coarseCorrFieldsPtr[cellI] *= sf;
                    }
                }
#ifdef SW_SLAVE
                else
                {
                    init_MVM_Arrays(&arrays1, nCells);
                    arrays1.A1Ptr = coarseCorrFieldsPtr;
                    arrays1.k1    = sf;
                    vectorOps_host(&arrays1, &slave_userFunc_aEk1Mua);
                }
#endif
            }

            //- correct the residual with the new solution
            agglomeration_.coarseMatrixLevels(leveli).spMV
            (
                ACf,
                coarseCorrFields[leveli]
            );

            scalar* coarseSourcesPtr = coarseSources[leveli].values();
            scalar* ACfPtr = ACf.values();
            IFNOT_SWACC
            {
                forAll(cellI, nCells)
                {
                    coarseSourcesPtr[cellI] -= ACfPtr[cellI];
                }
            }
#ifdef SW_SLAVE
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = coarseSourcesPtr;
                arrays1.A2Ptr = ACfPtr;
                vectorOps_host(&arrays1, &slave_userFunc_aEaMib);
            }
#endif
        }

        //- residual is equal to source
        //- note fineLevelIndex = leveli + 1, as leveli
        //  is coarse level and starts from 0
        agglomeration_.restrictField
        (
            coarseSources[leveli + 1],
            coarseSources[leveli],
            leveli + 1
        );
    }

    //- solve the coarsest level with either an iterative or direct solver

    if(coarsestLevel >= 0)
    {
        solveCoarsestLevel
        (
            coarseCorrFields[coarsestLevel],
            coarseSources[coarsestLevel]
        );
    }

    //- smoothing and prolongation of the coarse correction fields
    //  (going to finer levels)
    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        //- create a field for the pre-smoothed correction field
        //  as a sub-field of the finestCorrection which is not
        //  currently being used

        label nCells = coarseCorrFields[leveli].size();
        scalar* coarseCorrFieldsPtr = coarseCorrFields[leveli].values();

        scalarField preSmoothedCoarseCorrField(nCells);
        scalar* preSmoothedCoarseCorrFieldPtr = preSmoothedCoarseCorrField.values();

        //- only store the preSmoothedCoarseCorrField if pre-smoothing is used
        if (nPreSweeps_)
        {
            IFNOT_SWACC
            {
                forAll(cellI, nCells)
                {
                    preSmoothedCoarseCorrFieldPtr[cellI] = coarseCorrFieldsPtr[cellI];
                }
            }
#ifdef SW_SLAVE
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = preSmoothedCoarseCorrFieldPtr;
                arrays1.A2Ptr = coarseCorrFieldsPtr;
                vectorCopy_host(&arrays1);
            }
#endif
        }

        agglomeration_.prolongField
        (
            coarseCorrFields[leveli],
            coarseCorrFields[leveli + 1],
            leveli + 1
        );

        //- scale coarse-grid correction field
        //  but not on the coarsest level because it evaluates to 1
        if (scaleCorrection_ && leveli < coarsestLevel - 1)
        {
            //- create A.psi for this coarse level as a sub-field of Apsi
            scalarField ACf(nCells);

            scalar sf = scalingFactor
            (
                ACf,
                coarseCorrFields[leveli],
                coarseSources[leveli],
                agglomeration_.coarseMatrixLevels(leveli)
            );

            IFNOT_SWACC
            {
                forAll(cellI, nCells)
                {
                    coarseCorrFieldsPtr[cellI] *= sf;
                }
            }
#ifdef SW_SLAVE
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = coarseCorrFieldsPtr;
                arrays1.k1    = sf;
                vectorOps_host(&arrays1, &slave_userFunc_aEk1Mua);
            }
#endif
        }

        //- only add the preSmoothedCoarseCorrField if pre-smoothing is used
        if (nPreSweeps_)
        {
            IFNOT_SWACC
            {
                forAll(cellI, nCells)
                {
                    coarseCorrFieldsPtr[cellI] += preSmoothedCoarseCorrFieldPtr[cellI];
                }
            }
#ifdef SW_SLAVE
            else
            {
                init_MVM_Arrays(&arrays1, nCells);
                arrays1.A1Ptr = coarseCorrFieldsPtr;
                arrays1.A2Ptr = preSmoothedCoarseCorrFieldPtr;
                vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
            }
#endif
        }

        smoothers_[leveli+1].smooth
        (
            coarseCorrFields[leveli],
            agglomeration_.coarseMatrixLevels(leveli),
            coarseSources[leveli],
            nPostSweeps_ + leveli
        );
    }

    //- prolong the finest level correction

    if(coarsestLevel >= 0)
    {
        agglomeration_.prolongField
        (
            finestCorrection,
            coarseCorrFields[0],
            0
        );
    }

    if (scaleCorrection_)
    {
        //- calculate finest level scaling factor
        scalar fsf = scalingFactor
        (
            Apsi,
            finestCorrection,
            finestResidual,
            finestMatrix_
        );

        label nCells = psi.size();
        scalar* psiPtr = psi.values();
        scalar* finestCorrectionPtr = finestCorrection.values();
        IFNOT_SWACC
        {
            forAll(i, nCells)
            {
                psiPtr[i] += fsf*finestCorrectionPtr[i];
            }
        }
#ifdef SW_SLAVE
        else
        {
            MVM_Arrays arrays1;
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = psiPtr;
            arrays1.A2Ptr = finestCorrectionPtr;
            arrays1.k1    = fsf;
            // psi += fsf*finestCorrection;
            vectorOps_host(&arrays1, &slave_userFunc_aEaPk1Mub);
        }
#endif
    }
    else
    {
        label nCells = psi.size();
        scalar* psiPtr = psi.values();
        scalar* finestCorrectionPtr = finestCorrection.values();
        IFNOT_SWACC
        {
            forAll(i, nCells)
            {
                psiPtr[i] += finestCorrectionPtr[i];
            }
        }
#ifdef SW_SLAVE
        else
        {
            MVM_Arrays arrays1;
            init_MVM_Arrays(&arrays1, nCells);
            arrays1.A1Ptr = psiPtr;
            arrays1.A2Ptr = finestCorrectionPtr;
            // psi += finestCorrection;
            vectorOps_host(&arrays1, &slave_userFunc_aEaPb);
        }
#endif
    }

    smoothers_[0].smooth
    (
        psi,
        finestMatrix_,
        source,
        nFinestSweeps_
    );
}
