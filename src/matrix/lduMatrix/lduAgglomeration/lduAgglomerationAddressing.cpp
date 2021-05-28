#include <algorithm>
#include <vector>

#include "labelPair.hpp"
#include "lduAgglomeration.hpp"
#include "printUNAP.hpp"

void UNAP::lduAgglomeration::agglomerateLduAddressing(
    const label fineLevelIndex)
{
  const lduMatrix &fineA = matrixLevel(fineLevelIndex);

  const labelVector &upperAddr = fineA.upperAddr();
  const labelVector &lowerAddr = fineA.lowerAddr();

  label nFineFaces = upperAddr.size();

  //- get restriction map for current level
  const labelVector &restrictMap = restrictAddressing(fineLevelIndex);

#ifdef DEBUG
  if (restrictMap.size() != fineA.nCells())
  {
    commcator_->log()
        << "Error in agglomerateLduAddressing: restrict map does not "
           "correspond to fine level."
        << ENDL;
    commcator_->log() << "Sizes: restrictMap: " << restrictMap.size()
                      << " nEqns: " << fineA.nCells() << ENDL;
    ERROR_EXIT;
  }
#endif

  //- get the number of coarse cells
  const label nCoarseCells = nCells_[fineLevelIndex];

  //- storage for coarse cell neighbors and coefficients
  //- guess initial maximum number of neighbors in coarse cell
  label maxNnbrs = 10;

  //- number of faces for each coarse-cell
  labelVector cCellnFaces(nCoarseCells, 0, this->commcator_);

  //- setup initial packed storage for coarse-cell faces
  labelVector cCellFaces(maxNnbrs * nCoarseCells, this->commcator_);

  //- create face-restriction addressing
  labelVector *faceRestrictAddrPtr =
      new labelVector(nFineFaces, this->commcator_);
  labelVector &faceRestrictAddr = *faceRestrictAddrPtr;
  faceRestrictAddressing_.setLevel(fineLevelIndex, faceRestrictAddr);

  //- initial neighbor array (not in upper-triangle order)
  labelVector initCoarseNeighb(nFineFaces, this->commcator_);

  //- counter for coarse faces
  label nCoarseFaces = 0;

  //- loop through all fine faces
  forAll(fineFacei, upperAddr.size())
  {
    label rmUpperAddr = restrictMap[upperAddr[fineFacei]];
    label rmLowerAddr = restrictMap[lowerAddr[fineFacei]];

    if (rmUpperAddr == rmLowerAddr)
    {
      //- for each fine face inside of a coarse cell keep the address
      //  of the cell corresponding to the face in the faceRestrictAddr
      //  as a negative index
      faceRestrictAddr[fineFacei] = -(rmUpperAddr + 1);
    }
    else
    {
      //- this face is a part of a coarse face
      label cOwn = rmUpperAddr;
      label cNei = rmLowerAddr;

      //- get coarse owner and neighbor
      if (rmUpperAddr > rmLowerAddr)
      {
        cOwn = rmLowerAddr;
        cNei = rmUpperAddr;
      }

      //- check the neighbor to see if this face has already been found
      label *ccFaces = &cCellFaces[maxNnbrs * cOwn];

      bool nbrFound = false;
      label &ccnFaces = cCellnFaces[cOwn];

      for (int i = 0; i < ccnFaces; i++)
      {
        if (initCoarseNeighb[ccFaces[i]] == cNei)
        {
          nbrFound = true;
          faceRestrictAddr[fineFacei] = ccFaces[i];
          break;
        }
      }

      if (!nbrFound)
      {
        if (ccnFaces >= maxNnbrs)
        {
          label oldMaxNnbrs = maxNnbrs;
          maxNnbrs *= 2;

          cCellFaces.SET_size(maxNnbrs * nCoarseCells);

          for (label i = cCellnFaces.size() - 1; i >= 0; i--)
          {
            label *oldCcNbrs = &cCellFaces[oldMaxNnbrs * i];
            label *newCcNbrs = &cCellFaces[maxNnbrs * i];

            for (int j = 0; j < cCellnFaces[i]; j++)
            {
              newCcNbrs[j] = oldCcNbrs[j];
            }
          }

          ccFaces = &cCellFaces[maxNnbrs * cOwn];
        }

        ccFaces[ccnFaces] = nCoarseFaces;
        initCoarseNeighb[nCoarseFaces] = cNei;
        faceRestrictAddr[fineFacei] = nCoarseFaces;
        ccnFaces++;

        //- new coarse face created
        nCoarseFaces++;
      }
    }
  }  //- end for all fine faces

  //- renumber into upper-triangular order

  //- all coarse owner-neighbor storage
  labelVector coarseOwner(nCoarseFaces, commcator_);
  labelVector coarseNeighbour(nCoarseFaces, commcator_);
  labelVector coarseFaceMap(nCoarseFaces, commcator_);

  label coarseFacei = 0;

  forAll(cci, cCellnFaces.size())
  {
    label *cFaces = &cCellFaces[maxNnbrs * cci];
    label ccnFaces = cCellnFaces[cci];

    for (int i = 0; i < ccnFaces; i++)
    {
      coarseOwner[coarseFacei] = cci;
      coarseNeighbour[coarseFacei] = initCoarseNeighb[cFaces[i]];
      coarseFaceMap[cFaces[i]] = coarseFacei;
      coarseFacei++;
    }
  }

  forAll(fineFacei, faceRestrictAddr.size())
  {
    if (faceRestrictAddr[fineFacei] >= 0)
    {
      faceRestrictAddr[fineFacei] = coarseFaceMap[faceRestrictAddr[fineFacei]];
    }
  }

  //- clear the temporary storage for the coarse cell data
  cCellnFaces.SET_size(0);
  cCellFaces.SET_size(0);
  initCoarseNeighb.SET_size(0);
  coarseFaceMap.SET_size(0);

  //- Set the coarse ldu addressing onto the list
  lduMatrix *coarseAPtr = new lduMatrix(this->commcator_);
  lduMatrix &coarseA = *coarseAPtr;
  coarseMatrixLevels_.setLevel(fineLevelIndex, coarseA);

  coarseA.size(nCoarseCells);
  coarseA.SET_upperAddr(coarseNeighbour);
  coarseA.SET_lowerAddr(coarseOwner);

  if (this->commcator_->getMySize() > 1)
  {
    //- create coarse-level interfaces
    //- get reference to fine-level interfaces
    const interfaces &fineInterfaces = fineA.matrixInterfaces();
    const label interfacesSize = fineInterfaces.size();
    coarseA.matrixInterfaces(interfacesSize);
    interfaces &coarseInterfaces = coarseA.matrixInterfaces();

    //- create
    PtrList<labelVector> localAddressingList(interfacesSize);
    PtrList<labelVector> neighbourAddressingList(interfacesSize);

    string sendTaskName[interfacesSize];
    string recvTaskName[interfacesSize];
    //- initialize transfer of restrict addressing on the interface
    forAll(inti, interfacesSize)
    {
      const patch &finePatchI = fineInterfaces.patchList(inti);
      const label finePatchIFaceSize = finePatchI.size();

      //- point to PtrList, do not need to delete manually
      labelVector *localAddressingPtr =
          new labelVector(finePatchIFaceSize, this->commcator_);
      labelVector *neighbourAddressingPtr =
          new labelVector(finePatchIFaceSize, this->commcator_);

      labelVector &localAddressing = *localAddressingPtr;
      labelVector &neighbourAddressing = *neighbourAddressingPtr;

      localAddressingList.setLevel(inti, localAddressing);
      neighbourAddressingList.setLevel(inti, neighbourAddressing);

      const labelVector &faceCells = finePatchI.faceCells();
      forAll(faceI, finePatchIFaceSize)
      {
        localAddressing[faceI] = restrictMap[faceCells[faceI]];
      }

      const label neighbProcNo = finePatchI.neighbProcNo();

      char ch[128];

      sprintf(
          ch, "Send_%05d_Recv_%05d", this->commcator_->getMyId(), neighbProcNo);
      sendTaskName[inti] = ch;
      sprintf(
          ch, "Send_%05d_Recv_%05d", neighbProcNo, this->commcator_->getMyId());
      recvTaskName[inti] = ch;
      this->commcator_->send(sendTaskName[inti],
                             localAddressing.begin(),
                             finePatchIFaceSize * sizeof(label),
                             neighbProcNo);
      this->commcator_->recv(recvTaskName[inti],
                             neighbourAddressing.begin(),
                             finePatchIFaceSize * sizeof(label),
                             neighbProcNo);
    }

    forAll(inti, interfacesSize)
    {
      this->commcator_->finishTask(sendTaskName[inti]);
      this->commcator_->finishTask(recvTaskName[inti]);
      //- get coarse patch members: faceCells
      //- get patchFaceRestrictAddressing_

      const patch &finePatchI = fineInterfaces.patchList(inti);
      const label myProcNo = finePatchI.myProcNo();
      const label neighbProcNo = finePatchI.neighbProcNo();

      labelVector &localAddressing = localAddressingList[inti];
      labelVector &neighbourAddressing = neighbourAddressingList[inti];

      std::vector<labelPair> cellsToCoarseFace;
      std::vector<label> dynFaceCells;
      std::vector<label> dynFaceRestrictAddressing;

      forAll(ffi, localAddressing.size())
      {
        labelPair cellPair(-1, -1);

        //- do switching on master/slave indexes based on the owner/neighbor of
        //  the processor index such that both sides get the same answer.
        if (myProcNo < neighbProcNo)
        {
          //- master side
          cellPair.first(localAddressing[ffi]);
          cellPair.second(neighbourAddressing[ffi]);
        }
        else
        {
          //- slave side
          cellPair.first(neighbourAddressing[ffi]);
          cellPair.second(localAddressing[ffi]);
        }

        std::vector<labelPair>::iterator it = std::find(
            cellsToCoarseFace.begin(), cellsToCoarseFace.end(), cellPair);

        if (it == cellsToCoarseFace.end())
        {
          //- new coarse face
          label coarseI = dynFaceCells.size();
          dynFaceRestrictAddressing.push_back(coarseI);
          dynFaceCells.push_back(localAddressing[ffi]);
          cellPair.faceI(coarseI);
          cellsToCoarseFace.push_back(cellPair);
        }
        else
        {
          //- coarse face exists
          dynFaceRestrictAddressing.push_back(it->faceI());
        }
      }

      patch *coarsePatchIPtr =
          new patch(dynFaceCells.size(), myProcNo, neighbProcNo);
      patch &coarsePatchI = *coarsePatchIPtr;

      labelVector *coarseFaceCellsPtr =
          new labelVector(dynFaceCells.size(), this->commcator_);
      labelVector &coarseFaceCells = *coarseFaceCellsPtr;

      forAll(i, (label)dynFaceCells.size())
      {
        coarseFaceCells[i] = dynFaceCells[i];
      }

      labelVector *faceRestrictAddressingPtr =
          new labelVector(dynFaceRestrictAddressing.size(), this->commcator_);
      labelVector &faceRestrictAddressing = *faceRestrictAddressingPtr;

      forAll(i, (label)dynFaceRestrictAddressing.size())
      {
        faceRestrictAddressing[i] = dynFaceRestrictAddressing[i];
      }

      coarsePatchI.faceCells(coarseFaceCells);
      coarsePatchI.faceRestrictAddressing(faceRestrictAddressing);

      coarseInterfaces.patchList().setLevel(inti, coarsePatchI);
    }
  }
}
