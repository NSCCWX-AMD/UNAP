#include "lduAgglomeration.hpp"
#include "printUNAP.hpp"

#ifdef SW_SLAVE

void UNAP::lduAgglomeration::agglConstructRSSIterator()
{
  const label minFaces = SpMVAccSize;
  forAll(i, nCreatedLevels_)
  {
    const lduMatrix &fineM = matrixLevel(i);
    const label faceSize_f = fineM.upperAddr().size();

    if (faceSize_f >= minFaces)
    {
      lduMatrix &coarseM = coarseMatrixLevels_[i];
      const label faceSize_c = coarseM.upperAddr().size();
      if (faceSize_c >= minFaces)
      {
        coarseM.constructRSSIterator();
      }
    }
  }
}

void UNAP::lduAgglomeration::agglomerationReorderTopo()
{
  const label minFaces = SpMVAccSize;
  forAll(i, nCreatedLevels_)
  {
    const lduMatrix &fineM = matrixLevel(i);
    const label cellSize_f = fineM.size();
    const label faceSize_f = fineM.upperAddr().size();

    if (faceSize_f >= minFaces)
    {
      const label *unatCellMap_f = fineM.unatCellMap();
      const label *unatFaceMap_f = fineM.unatEdgeMap();

      //- cell
      labelVector resMapOld(restrictAddressing_[i]);
      label *resMapPtr_new = restrictAddressing_[i].begin();
      label *resMapPtr_old = resMapOld.begin();

      //- face
      labelVector faceResMapOld(faceRestrictAddressing_[i]);
      label *faceResMapPtr_new = faceRestrictAddressing_[i].begin();
      label *faceResMapPtr_old = faceResMapOld.begin();

      lduMatrix &coarseM = coarseMatrixLevels_[i];
      const label faceSize_c = coarseM.upperAddr().size();

      if (faceSize_c >= minFaces)
      {
        coarseM.constructMLBIterator();

        //- cell
        const label *unatCellMap_c = coarseM.unatCellMap();

        forAll(j, cellSize_f)
        {
          resMapPtr_new[unatCellMap_f[j]] = unatCellMap_c[resMapPtr_old[j]];
        }

        //- face
        const label *unatFaceMap_c = coarseM.unatEdgeMap();

        forAll(j, faceSize_f)
        {
          if (faceResMapPtr_old[j] >= 0)
          {
            faceResMapPtr_new[abs(unatFaceMap_f[j]) - 1] =
                abs(unatFaceMap_c[faceResMapPtr_old[j]]) - 1;
          }
          else
          {
            label posOld = -1 - faceResMapPtr_old[j];
            faceResMapPtr_new[abs(unatFaceMap_f[j]) - 1] =
                -(unatCellMap_c[posOld] + 1);
          }
        }
      }
      else
      {
        //- cell
        forAll(j, cellSize_f)
        {
          resMapPtr_new[unatCellMap_f[j]] = resMapPtr_old[j];
        }

        //- face
        forAll(j, faceSize_f)
        {
          faceResMapPtr_new[abs(unatFaceMap_f[j]) - 1] = faceResMapPtr_old[j];
        }
      }
    }
  }
}

#endif
