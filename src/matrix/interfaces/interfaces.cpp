#include "interfaces.hpp"

UNAP::interfaces::interfaces(PtrList<patch> &patches, Communicator *other_comm)
    : patches_(patches),
      sendBuffer_(NULL),
      recvBuffer_(NULL),
      sendTaskName_(NULL),
      recvTaskName_(NULL),
      locPosition_(NULL),
      destRank_(NULL),
      commcator_(other_comm)
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
  if (!commcator_)
    commcator_ = psi.getCommunicator();
  else if (commcator_ != psi.getCommunicator())
  {
    commcator_->log()
        << "Error" << __FILE__ << " " << __LINE__
        << "The communicators between interfaces and Apsi are different\n";
    ERROR_EXIT;
  }

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
    forAll(faceI, locSize)
    {
      sendBuffer_[faceI + locPosition_[i]] = psiPtr[faceCellsPtr[faceI]];
    }

    char ch[128];

    sprintf(
        ch, "Send_%05d_Recv_%05d", this->commcator_->getMyId(), destRank_[i]);
    sendTaskName_[i] = ch;
    sprintf(
        ch, "Send_%05d_Recv_%05d", destRank_[i], this->commcator_->getMyId());
    recvTaskName_[i] = ch;

    commcator_->send(sendTaskName_[i],
                     &(sendBuffer_[locPosition_[i]]),
                     sizeof(scalar) * locSize,
                     destRank_[i]);
    commcator_->recv(recvTaskName_[i],
                     &(recvBuffer_[locPosition_[i]]),
                     sizeof(scalar) * locSize,
                     destRank_[i]);
  }
}

void UNAP::interfaces::updateMatrixInterfaces(scalarVector &Apsi) const
{
  if (commcator_ != Apsi.getCommunicator())
  {
    commcator_->log()
        << "Error" << __FILE__ << " " << __LINE__
        << "The communicators between interfaces and Apsi are different\n";
    ERROR_EXIT;
  }
  scalar *ApsiPtr = Apsi.begin();
  label numInterfaces = patches_.size();

  forAll(i, numInterfaces)
  {
    commcator_->finishTask(sendTaskName_[i]);
    commcator_->finishTask(recvTaskName_[i]);

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
