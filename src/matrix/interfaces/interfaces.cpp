#include "interfaces.hpp"

UNAP::interfaces::interfaces(PtrList<patch>& patches)
:
	patches_(patches),
	sendBuffer_(NULL),
	recvBuffer_(NULL),
	sendRecvRequests_(NULL),
	locPosition_(NULL),
	destRank_(NULL)
{}


UNAP::interfaces::~interfaces()
{
	delete &patches_;
	DELETE_POINTER(sendBuffer_);
	DELETE_POINTER(recvBuffer_);
	DELETE_POINTER(sendRecvRequests_);
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
	sendRecvRequests_ = new MPI_Request[2*numInterfaces];

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

		MPI_Isend
		(
			&(sendBuffer_[locPosition_[i]]),
			locSize,
			MPI_SCALAR,
			destRank_[i],
			1,
			MPI_COMM_WORLD,
			&sendRecvRequests_[i]
		);

		MPI_Irecv
		(
			&(recvBuffer_[locPosition_[i]]),
			locSize,
			MPI_SCALAR,
			destRank_[i],
			1,
			MPI_COMM_WORLD,
			&sendRecvRequests_[i + numInterfaces]
		);
	}
}


void UNAP::interfaces::updateMatrixInterfaces
(
    scalarField& Apsi
) const
{
	scalar* ApsiPtr = Apsi.begin();
	label numInterfaces = patches_.size();
	MPI_Waitall(2*numInterfaces, &sendRecvRequests_[0], MPI_STATUSES_IGNORE);

	forAll(i, numInterfaces)
	{
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
	DELETE_POINTER(sendRecvRequests_)
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
