#include "lduMatrix.hpp"

#ifdef SW_SLAVE
#include "SW_UNAPInterface.h"
#endif

#define PARRUN (commcator_->getMySize() > 1)
#define MYID (commcator_->getMyId())

UNAP::lduMatrix::lduMatrix(Communicator *other_comm)
    : nCells_(-1),
      lowerAddrPtr_(NULL),
      upperAddrPtr_(NULL),
      lowerPtr_(NULL),
      diagPtr_(NULL),
      upperPtr_(NULL),
      interfacesPtr_(NULL),
      losortPtr_(NULL),
      ownerStartPtr_(NULL),
      losortStartPtr_(NULL),
      matrix(other_comm)
{
#ifdef SW_SLAVE
  mlbIter_ = NULL;
  rssIter_ = NULL;
  unatIter_ = NULL;
#endif
}

UNAP::lduMatrix::lduMatrix(const label &nCells,
                           const labelVector &lowerAddr,
                           const labelVector &upperAddr,
                           const scalarVector &lower,
                           const scalarVector &diag,
                           const scalarVector &upper,
                           Communicator *other_comm)
    : nCells_(nCells),
      lowerAddrPtr_(NULL),
      upperAddrPtr_(NULL),
      lowerPtr_(NULL),
      diagPtr_(NULL),
      upperPtr_(NULL),
      interfacesPtr_(NULL),
      losortPtr_(NULL),
      ownerStartPtr_(NULL),
      losortStartPtr_(NULL),
      matrix(other_comm)
{
  ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelVector)
  ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelVector)
  ALLOCATE_POINTER(diagPtr_, diag, scalarVector)
  ALLOCATE_POINTER(upperPtr_, upper, scalarVector)

  if ((&lower) == (&upper))
  {
    lowerPtr_ = upperPtr_;
  }
  else
  {
    ALLOCATE_POINTER(lowerPtr_, lower, scalarVector)
  }

#ifdef SW_SLAVE
  mlbIter_ = NULL;
  rssIter_ = NULL;
  unatIter_ = NULL;
#endif
}

UNAP::lduMatrix::lduMatrix(const label &nCells,
                           const labelVector &lowerAddr,
                           const labelVector &upperAddr,
                           const scalarVector &lower,
                           const scalarVector &diag,
                           const scalarVector &upper,
                           const bool reUse,
                           Communicator *other_comm)
    : nCells_(nCells),
      lowerAddrPtr_(NULL),
      upperAddrPtr_(NULL),
      lowerPtr_(NULL),
      diagPtr_(NULL),
      upperPtr_(NULL),
      interfacesPtr_(NULL),
      losortPtr_(NULL),
      ownerStartPtr_(NULL),
      losortStartPtr_(NULL),
      matrix(other_comm)
{
  if (reUse)
  {
    lowerAddrPtr_ = const_cast<labelVector *>(&lowerAddr);
    upperAddrPtr_ = const_cast<labelVector *>(&upperAddr);
    diagPtr_ = const_cast<scalarVector *>(&diag);
    lowerPtr_ = const_cast<scalarVector *>(&lower);
    upperPtr_ = const_cast<scalarVector *>(&upper);
  }
  else
  {
    ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelVector)
    ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelVector)
    ALLOCATE_POINTER(diagPtr_, diag, scalarVector)
    ALLOCATE_POINTER(upperPtr_, upper, scalarVector)

    if ((&lower) == (&upper))
    {
      lowerPtr_ = upperPtr_;
    }
    else
    {
      ALLOCATE_POINTER(lowerPtr_, lower, scalarVector)
    }
  }

#ifdef SW_SLAVE
  mlbIter_ = NULL;
  rssIter_ = NULL;
  unatIter_ = NULL;
#endif
}

UNAP::lduMatrix::lduMatrix(const label &nCells,
                           const labelVector &lowerAddr,
                           const labelVector &upperAddr,
                           Communicator *other_comm)
    : nCells_(nCells),
      lowerAddrPtr_(NULL),
      upperAddrPtr_(NULL),
      lowerPtr_(NULL),
      diagPtr_(NULL),
      upperPtr_(NULL),
      interfacesPtr_(NULL),
      losortPtr_(NULL),
      ownerStartPtr_(NULL),
      losortStartPtr_(NULL),
      matrix(other_comm)
{
  ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelVector)
  ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelVector)

#ifdef SW_SLAVE
  mlbIter_ = NULL;
  rssIter_ = NULL;
  unatIter_ = NULL;
#endif
}

UNAP::lduMatrix::lduMatrix(const label &nCells,
                           const labelVector &lowerAddr,
                           const labelVector &upperAddr,
                           const bool reUse,
                           Communicator *other_comm)
    : nCells_(nCells),
      lowerAddrPtr_(NULL),
      upperAddrPtr_(NULL),
      lowerPtr_(NULL),
      diagPtr_(NULL),
      upperPtr_(NULL),
      interfacesPtr_(NULL),
      losortPtr_(NULL),
      ownerStartPtr_(NULL),
      losortStartPtr_(NULL),
      matrix(other_comm)
{
  if (reUse)
  {
    lowerAddrPtr_ = const_cast<labelVector *>(&lowerAddr);
    upperAddrPtr_ = const_cast<labelVector *>(&upperAddr);
  }
  else
  {
    ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelVector)
    ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelVector)
  }

#ifdef SW_SLAVE
  mlbIter_ = NULL;
  rssIter_ = NULL;
  unatIter_ = NULL;
#endif
}

UNAP::lduMatrix::~lduMatrix()
{
  if (lowerPtr_ != upperPtr_)
  {
    DELETE_OBJECT_POINTER(upperPtr_)
    DELETE_OBJECT_POINTER(lowerPtr_)
  }
  else
  {
    DELETE_OBJECT_POINTER(upperPtr_)
    lowerPtr_ = NULL;
  }

  DELETE_OBJECT_POINTER(lowerAddrPtr_)
  DELETE_OBJECT_POINTER(upperAddrPtr_)
  DELETE_OBJECT_POINTER(diagPtr_)
  DELETE_OBJECT_POINTER(interfacesPtr_)
  DELETE_OBJECT_POINTER(losortPtr_)
  DELETE_OBJECT_POINTER(ownerStartPtr_)
  DELETE_OBJECT_POINTER(losortStartPtr_)

#ifdef SW_SLAVE
  DELETE_OBJECT_POINTER(mlbIter_)
  DELETE_OBJECT_POINTER(rssIter_)
  unatIter_ = NULL;
#endif
}

UNAP::labelVector &UNAP::lduMatrix::lowerAddr() const
{
#ifdef DEBUG
  CHECK_POINTER(lowerAddrPtr_)
#endif

  return *lowerAddrPtr_;
}

UNAP::labelVector &UNAP::lduMatrix::upperAddr() const
{
#ifdef DEBUG
  CHECK_POINTER(upperAddrPtr_)
#endif

  return *upperAddrPtr_;
}

UNAP::scalarVector &UNAP::lduMatrix::lower() const
{
#ifdef DEBUG
  CHECK_POINTER(lowerPtr_)
#endif

  return *lowerPtr_;
}

UNAP::scalarVector &UNAP::lduMatrix::diag() const
{
#ifdef DEBUG
  CHECK_POINTER(diagPtr_)
#endif

  return *diagPtr_;
}

UNAP::scalarVector &UNAP::lduMatrix::upper() const
{
#ifdef DEBUG
  CHECK_POINTER(upperPtr_)
#endif

  return *upperPtr_;
}

void UNAP::lduMatrix::calcOwnerStart() const
{
  const labelVector &own = lowerAddr();

  ownerStartPtr_ = new labelVector(nCells() + 1, own.size(), this->commcator_);

  labelVector &ownStart = *ownerStartPtr_;

  //- set up first lookup by hand
  ownStart[0] = 0;
  label nOwnStart = 0;
  label i = 1;

  forAll(faceI, own.size())
  {
    label curOwn = own[faceI];

    if (curOwn > nOwnStart)
    {
      while (i <= curOwn)
      {
        ownStart[i++] = faceI;
      }

      nOwnStart = curOwn;
    }
  }
}

const UNAP::labelVector &UNAP::lduMatrix::ownerStartAddr() const
{
  if (!ownerStartPtr_)
  {
    calcOwnerStart();
  }

  return *ownerStartPtr_;
}

void UNAP::lduMatrix::calcLosort() const
{
  //- scan the neighbor list to find out how many times the cell
  //  appears as a neighbor of the face. Done this way to avoid guessing
  //  and resizing list
  labelVector nNbrOfFace(nCells(), 0);

  const labelVector &nbr = upperAddr();

  forAll(nbrI, nbr.size()) { nNbrOfFace[nbr[nbrI]]++; }

  //- create temporary neighbor addressing
  label **cellNbrFaces = new label *[nCells_];
  forAll(cellI, nCells())
  {
    cellNbrFaces[cellI] = new label[nNbrOfFace[cellI]];
  }

  //- reset the list of number of neighbors to zero
  nNbrOfFace = 0;

  //- scatter the neighbor faces
  forAll(nbrI, nbr.size())
  {
    cellNbrFaces[nbr[nbrI]][nNbrOfFace[nbr[nbrI]]] = nbrI;

    nNbrOfFace[nbr[nbrI]]++;
  }

  //- gather the neighbors into the losort array
  losortPtr_ = new labelVector(nbr.size(), -1, this->commcator_);

  labelVector &lst = *losortPtr_;

  //- set counter for losort
  label lstI = 0;

  forAll(cellI, nCells_)
  {
    forAll(curNbrI, nNbrOfFace[cellI])
    {
      lst[lstI] = cellNbrFaces[cellI][curNbrI];
      lstI++;
    }
  }

  //- free
  forAll(cellI, nCells_)
  {
    delete[] cellNbrFaces[cellI];
    cellNbrFaces[cellI] = NULL;
  }
  delete[] cellNbrFaces;
  cellNbrFaces = NULL;
}

const UNAP::labelVector &UNAP::lduMatrix::losortAddr() const
{
  if (!losortPtr_)
  {
    calcLosort();
  }

  return *losortPtr_;
}

void UNAP::lduMatrix::calcLosortStart() const
{
  losortStartPtr_ = new labelVector(nCells() + 1, 0, this->commcator_);

  labelVector &lsrtStart = *losortStartPtr_;

  const labelVector &nbr = upperAddr();

  const labelVector &lsrt = losortAddr();

  //- set up first lookup by hand
  lsrtStart[0] = 0;
  label nLsrtStart = 0;
  label i = 0;

  forAll(faceI, lsrt.size())
  {
    //- get neighbor
    const label curNbr = nbr[lsrt[faceI]];

    if (curNbr > nLsrtStart)
    {
      while (i <= curNbr)
      {
        lsrtStart[i++] = faceI;
      }

      nLsrtStart = curNbr;
    }
  }

  //- set up last lookup by hand
  lsrtStart[nCells()] = nbr.size();
}

const UNAP::labelVector &UNAP::lduMatrix::losortStartAddr() const
{
  if (!losortStartPtr_)
  {
    calcLosortStart();
  }

  return *losortStartPtr_;
}

void UNAP::lduMatrix::createInterfacesTopology(const label nNeiProcs,
                                               const label *destRank,
                                               const label *offDiagRows,
                                               const label *offDiagStarts)
{
  PtrList<patch> *patchesPtr = new PtrList<patch>(nNeiProcs);

  forAll(intI, nNeiProcs)
  {
    const label neighborID = destRank[intI] - 1;
    const label localSize = offDiagStarts[intI + 1] - offDiagStarts[intI];
    patch *patchIPtr = new patch(localSize, MYID, neighborID);

    label *localFaceCells = new label[localSize];

    forAll(faceI, localSize)
    {
      label start = offDiagStarts[intI] + faceI - 1;
      localFaceCells[faceI] = offDiagRows[start] - 1;
    }

    labelVector *locFaceCellsPtr =
        new labelVector(localFaceCells, localSize, true, this->commcator_);
    scalarVector *localDataPtr = new scalarVector(localSize, this->commcator_);

    patchIPtr->faceCells(*locFaceCellsPtr);
    patchIPtr->patchCoeffs(*localDataPtr);

    patchesPtr->setLevel(intI, *patchIPtr);
  }

  interfaces *interfacesLocalPtr =
      new interfaces(*patchesPtr, this->commcator_);
  matrixInterfaces(*interfacesLocalPtr);
}

void UNAP::lduMatrix::fillInterfacesCofficients(const label *offDiagStarts,
                                                const scalar *offDiagCoeffs)
{
  interfaces &lduAInter = matrixInterfaces();

  const label nNeiProcs = lduAInter.size();

  forAll(intI, nNeiProcs)
  {
    patch &patchI = lduAInter.patchList(intI);
    const label localSize = offDiagStarts[intI + 1] - offDiagStarts[intI];
    scalar *localData = patchI.patchCoeffs().begin();

    forAll(faceI, localSize)
    {
      label start = offDiagStarts[intI] + faceI - 1;
      localData[faceI] = offDiagCoeffs[start];
    }
  }
}

void UNAP::lduMatrix::initInterfaces(const scalarVector &psi) const
{
  if (PARRUN)
  {
    matrixInterfaces().initMatrixInterfaces(psi);
  }
}

void UNAP::lduMatrix::updateInterfaces(scalarVector &Apsi) const
{
  if (PARRUN)
  {
    matrixInterfaces().updateMatrixInterfaces(Apsi);
  }
}

void UNAP::lduMatrix::setMatrixCoeffients(const scalarVector &diag,
                                          const scalarVector &upper,
                                          const scalarVector &lower)
{
  nCells_ = diag.size();
  ALLOCATE_POINTER(diagPtr_, diag, scalarVector)

  ALLOCATE_POINTER(upperPtr_, upper, scalarVector)

  if ((&lower) == (&upper))
  {
    lowerPtr_ = upperPtr_;
  }
  else
  {
    ALLOCATE_POINTER(lowerPtr_, lower, scalarVector)
  }
}

void UNAP::lduMatrix::setMatrixCoeffients(const scalarVector &diag,
                                          const scalarVector &upper,
                                          const scalarVector &lower,
                                          const bool reUse)
{
  nCells_ = diag.size();
  if (reUse)
  {
    DELETE_OBJECT_POINTER(diagPtr_)
    this->diagPtr_ = const_cast<scalarVector *>(&diag);

    DELETE_OBJECT_POINTER(upperPtr_)
    this->upperPtr_ = const_cast<scalarVector *>(&upper);

    if ((&lower) == (&upper))
    {
      this->lowerPtr_ = this->upperPtr_;
    }
    else
    {
      DELETE_OBJECT_POINTER(lowerPtr_)
      this->lowerPtr_ = const_cast<scalarVector *>(&lower);
    }
  }
  else
  {
    setMatrixCoeffients(diag, upper, lower);
  }
}

void UNAP::lduMatrix::setMatrixTopology(const labelVector &upperAddr,
                                        const labelVector &lowerAddr)
{
  ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelVector)

  ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelVector)
}

void UNAP::lduMatrix::setMatrixTopology(const labelVector &upperAddr,
                                        const labelVector &lowerAddr,
                                        const bool reUse)
{
  if (reUse)
  {
    DELETE_OBJECT_POINTER(upperAddrPtr_)
    this->upperAddrPtr_ = const_cast<labelVector *>(&upperAddr);

    DELETE_OBJECT_POINTER(lowerAddrPtr_)
    this->lowerAddrPtr_ = const_cast<labelVector *>(&lowerAddr);
  }
  else
  {
    setMatrixTopology(upperAddr, lowerAddr);
  }
}

#ifdef SW_SLAVE

void UNAP::lduMatrix::constructMLBIterator()
{
  label *uPtr = this->upperAddr().begin();
  label *lPtr = this->lowerAddr().begin();
  label nFaces = this->lowerAddr().size();
  label nCells = this->size();

  if (nFaces >= SpMVAccSize)
  {
    label *cellWeights = new label[nCells];
    label *faceWeights = new label[nFaces];

    for (int i = 0; i < nCells; i++)
    {
      cellWeights[i] = 3;
    }
    for (int i = 0; i < nFaces; i++)
    {
      faceWeights[i] = 2;
    }
    mlbIter_ = constructMLBIteratorFromUNAP(
        lPtr, uPtr, cellWeights, faceWeights, nFaces);

    this->unatEdgeMap_ = getEdgeMap(this->mlbIter_);
    this->unatCellMap_ = getCellMap(this->mlbIter_);

    this->mlbIter_->reorderEdges(lPtr, uPtr, nFaces, nCells);
    if (PARRUN)
    {
      this->interfacesPtr_->reorderIntFaceCells(unatCellMap_);
    }
    DELETE_POINTER(cellWeights)
    DELETE_POINTER(faceWeights)
  }

  unatIter_ = mlbIter_;
}

void UNAP::lduMatrix::reorderLDUValues()
{
  scalar *uPtr = this->upper().begin();
  scalar *lPtr = this->lower().begin();
  scalar *dPtr = this->diag().begin();
  label nFaces = this->lower().size();
  label nCells = this->size();

  if (mlbIter_)
  {
    Arrays selfConnData, frontEdgeData, backEdgeData;
    constructSingleArray(selfConnData, 1, nCells, COPYIN, (swFloat *)dPtr);
    constructSingleArray(frontEdgeData, 1, nFaces, COPYIN, (swFloat *)uPtr);
    constructSingleArray(backEdgeData, 1, nFaces, COPYIN, (swFloat *)lPtr);
    this->mlbIter_->reorderVertexData(&selfConnData);
    this->mlbIter_->reorderEdgeData(&backEdgeData, &frontEdgeData);
  }
}

void UNAP::lduMatrix::restoreVector(scalarVector &vv)
{
  scalar *vvPtr = vv.begin();
  label nCells = vv.size();

  if (mlbIter_)
  {
    Arrays selfConnData;
    constructSingleArray(selfConnData, 1, nCells, COPYIN, (swFloat *)vvPtr);
    this->mlbIter_->restoreVertexData(&selfConnData);
  }
}

void UNAP::lduMatrix::reorderVector(scalarVector &vv)
{
  scalar *vvPtr = vv.begin();
  label nCells = vv.size();

  if (mlbIter_)
  {
    Arrays selfConnData;
    constructSingleArray(selfConnData, 1, nCells, COPYIN, (swFloat *)vvPtr);
    this->mlbIter_->reorderVertexData(&selfConnData);
  }
}

void UNAP::lduMatrix::constructRSSIterator()
{
  label *uPtr = this->upperAddr().begin();
  label *lPtr = this->lowerAddr().begin();
  label nFaces = this->lowerAddr().size();
  label nCells = this->size();

  if (nFaces >= SpMVAccSize)
  {
    label *cellWeights = new label[nCells];
    label *faceWeights = new label[nFaces];

    for (int i = 0; i < nCells; i++)
    {
      cellWeights[i] = 3;
    }
    for (int i = 0; i < nFaces; i++)
    {
      faceWeights[i] = 2;
    }
    rssIter_ = constructRSSIteratorFromUNAP(
        lPtr, uPtr, cellWeights, faceWeights, nFaces);

    DELETE_POINTER(cellWeights)
    DELETE_POINTER(faceWeights)
  }

  unatIter_ = rssIter_;
}

void UNAP::lduMatrix::setMlbIter(UNAT::MultiLevelBlockIterator *mlbIter)
{
  mlbIter_ = mlbIter;
  label fsize = lowerAddrPtr_->size();
  unatEdgeMap_ = new label[fsize];
  forAll(ii, fsize) { unatEdgeMap_[ii] = ii + 1; }

  unatCellMap_ = new label[nCells_];
  forAll(ii, nCells_) { unatCellMap_[ii] = ii; }
}

#endif
