#include "interfaces.hpp"

UNAP::interfaces::interfaces(PtrList<patch>& patches)
:
	patches_(patches),
	sendBuffer_(NULL),
	recvBuffer_(NULL),
	sendRecvTaskName_(NULL),
	locPosition_(NULL),
	destRank_(NULL)
{}


UNAP::interfaces::~interfaces()
{
	delete &patches_;
	DELETE_POINTER(sendBuffer_);
	DELETE_POINTER(recvBuffer_);
	DELETE_POINTER(sendRecvTaskName_);
	DELETE_POINTER(locPosition_);
	DELETE_POINTER(destRank_);
}


void UNAP::interfaces::initMatrixInterfaces
(
    const scalarField& psi
) const
{
	label numInterfaces = patches_.size();
	locPosition_ = new label[numInterfaces+1];
	destRank_ = new label[numInterfaces];
	sendRecvTaskName_ = new string[numInterfaces];

	locPosition_[0] = 0;
	forAll(i, numInterfaces)
	{
		patch& patchI = patches_[i];
		destRank_[i] = patchI.neighbProcNo();

		locPosition_[i+1] = patchI.size() + locPosition_[i];
	}

	label bufferSize = locPosition_[numInterfaces];

	sendBuffer_ = new scalar[bufferSize];
	recvBuffer_ = new scalar[bufferSize];


	const scalar* const psiPtr = psi.begin();

	forAll(i, numInterfaces)
	{
		patch& patchI = patches_[i];
		const label* const faceCellsPtr = patchI.faceCells().begin();

		label locSize = locPosition_[i+1] - locPosition_[i];
		forAll(faceI, locSize)
		{
			sendBuffer_[faceI+locPosition_[i]] = psiPtr[faceCellsPtr[faceI]];
		}

		char ch[8];
		sendRecvTaskName_[i] = "sendRecv_";
		sprintf(ch,"%05%d",i);
		sendRecvTaskName_[i] += ch;

		UNAP::unapMPI::unapCommunicator().send(sendRecvTaskName_[i], &(sendBuffer_[locPosition_[i]]), sizeof(scalar)*locSize,destRank_[i]);
		UNAP::unapMPI::unapCommunicator().recv(sendRecvTaskName_[i], &(recvBuffer_[locPosition_[i]]), sizeof(scalar)*locSize,destRank_[i]);
	}
}


void UNAP::interfaces::updateMatrixInterfaces
(
    scalarField& Apsi
) const
{
	scalar* ApsiPtr = Apsi.begin();
	label numInterfaces = patches_.size();

	forAll(i, numInterfaces)
	{
		UNAP::unapMPI::unapCommunicator().finishTask(sendRecvTaskName_[i]);

		patch& patchI = patches_[i];
		const label* const faceCellsPtr = patchI.faceCells().begin();
		label locSize = locPosition_[i+1] - locPosition_[i];

		forAll(faceI, locSize)
		{
			ApsiPtr[faceCellsPtr[faceI]] +=
	        	patchI.patchCoeffs(faceI) * recvBuffer_[faceI+locPosition_[i]];
		}
	}

	DELETE_POINTER(sendBuffer_)
	DELETE_POINTER(recvBuffer_)
	DELETE_POINTER(sendRecvTaskName_)
	DELETE_POINTER(locPosition_)
	DELETE_POINTER(destRank_)
}


void UNAP::interfaces::reorderIntFaceCells(const label* cellMap)
{
	label numInterfaces = patches_.size();
	forAll(i, numInterfaces)
	{
		patch& patchI = patches_[i];
		patchI.reorderPatchFaceCells(cellMap);
	}
}
