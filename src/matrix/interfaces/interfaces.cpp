#include "interfaces.hpp"

UNAP::interfaces::interfaces(PtrList<patch> &patches)
    : patches_(patches),
      sendBuffer_(NULL),
      recvBuffer_(NULL),
      sendTaskName_(NULL),
      recvTaskName_(NULL),
      locPosition_(NULL),
      destRank_(NULL)
{
}

UNAP::interfaces::~interfaces()
{
  delete &patches_;
  DELETE_POINTER(sendBuffer_);
  DELETE_POINTER(recvBuffer_);
  DELETE_POINTER(sendTaskName_);
  DELETE_POINTER(recvTaskName_);
  DELETE_POINTER(locPosition_);
  DELETE_POINTER(destRank_);
}

void UNAP::interfaces::initMatrixInterfaces(const scalarVector &psi) const
{
  label numInterfaces = patches_.size();
  locPosition_ = new label[numInterfaces + 1];
  destRank_ = new label[numInterfaces];
  sendTaskName_ = new string[numInterfaces];
  recvTaskName_ = new string[numInterfaces];
  locPosition_[0] = 0;
  forAll(i, numInterfaces)
  {
    patch &patchI = patches_[i];
    destRank_[i] = patchI.neighbProcNo();

    locPosition_[i + 1] = patchI.size() + locPosition_[i];
  }

  label bufferSize = locPosition_[numInterfaces];

  sendBuffer_ = new scalar[bufferSize];
  recvBuffer_ = new scalar[bufferSize];

  const scalar *const psiPtr = psi.begin();

  forAll(i, numInterfaces)
  {
    patch &patchI = patches_[i];
    const label *const faceCellsPtr = patchI.faceCells().begin();

    label locSize = locPosition_[i + 1] - locPosition_[i];

    if (MYID == destRank_[i])
    {
      const label *const faceCells2Ptr = patchI.faceCells2().begin();
      forAll(faceI, locSize)
      {
        sendBuffer_[faceI + locPosition_[i]] = psiPtr[faceCells2Ptr[faceI]];
      }
    }
    else
    {
      forAll(faceI, locSize)
      {
        sendBuffer_[faceI + locPosition_[i]] = psiPtr[faceCellsPtr[faceI]];
      }
    }

    char ch[128];
    sprintf(ch, "Send_%05d_Recv_%05d", MYID, destRank_[i]);
    sendTaskName_[i] = ch;
    sprintf(ch, "Send_%05d_Recv_%05d", destRank_[i], MYID);
    recvTaskName_[i] = ch;

    UNAP::unapMPI::unapCommunicator().send(sendTaskName_[i],
                                           &(sendBuffer_[locPosition_[i]]),
                                           sizeof(scalar) * locSize,
                                           destRank_[i]);
    UNAP::unapMPI::unapCommunicator().recv(recvTaskName_[i],
                                           &(recvBuffer_[locPosition_[i]]),
                                           sizeof(scalar) * locSize,
                                           destRank_[i]);
  }
}

void UNAP::interfaces::updateMatrixInterfaces(scalarVector &Apsi) const
{
  scalar *ApsiPtr = Apsi.begin();
  label numInterfaces = patches_.size();

  forAll(i, numInterfaces)
  {
    UNAP::unapMPI::unapCommunicator().finishTask(sendTaskName_[i]);
    UNAP::unapMPI::unapCommunicator().finishTask(recvTaskName_[i]);

    patch &patchI = patches_[i];
    const label *const faceCellsPtr = patchI.faceCells().begin();
    label locSize = locPosition_[i + 1] - locPosition_[i];

    forAll(faceI, locSize)
    {
      ApsiPtr[faceCellsPtr[faceI]] +=
          patchI.patchCoeffs(faceI) * recvBuffer_[faceI + locPosition_[i]];
    }
  }

  DELETE_POINTER(sendBuffer_)
  DELETE_POINTER(recvBuffer_)
  DELETE_POINTER(sendTaskName_)
  DELETE_POINTER(recvTaskName_)
  DELETE_POINTER(locPosition_)
  DELETE_POINTER(destRank_)
}

void UNAP::interfaces::reorderIntFaceCells(const label *cellMap)
{
  label numInterfaces = patches_.size();
  forAll(i, numInterfaces)
  {
    patch &patchI = patches_[i];
    patchI.reorderPatchFaceCells(cellMap);
  }
}
