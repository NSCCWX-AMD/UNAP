#include "lduMatrix.hpp"

#ifdef SW_SLAVE
#include "SW_UNAPInterface.h"
#endif

UNAP::lduMatrix::lduMatrix()
:
    nCells_(-1),
    lowerAddrPtr_(NULL),
    upperAddrPtr_(NULL),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    interfacesPtr_(NULL),
    losortPtr_(NULL),
    ownerStartPtr_(NULL),
    losortStartPtr_(NULL)
{
#ifdef SW_SLAVE
    mlbIter_ = NULL;
    rssIter_ = NULL;
    unatIter_ = NULL;
#endif
}


UNAP::lduMatrix::lduMatrix
(
	const label&       nCells,
    const labelField&  lowerAddr,
    const labelField&  upperAddr,
    const scalarField& lower,
    const scalarField& diag,
    const scalarField& upper
)
:
	nCells_(nCells),
	lowerAddrPtr_(NULL),
    upperAddrPtr_(NULL),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    interfacesPtr_(NULL),
    losortPtr_(NULL),
    ownerStartPtr_(NULL),
    losortStartPtr_(NULL)
{
    ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelField)
    ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelField)
    ALLOCATE_POINTER(diagPtr_,      diag,      scalarField)
    ALLOCATE_POINTER(upperPtr_,     upper,     scalarField)

    if((&lower) == (&upper))
    {
        lowerPtr_ = upperPtr_;
    }
    else
    {
        ALLOCATE_POINTER(lowerPtr_, lower, scalarField)
    }

#ifdef SW_SLAVE
    mlbIter_ = NULL;
    rssIter_ = NULL;
    unatIter_ = NULL;
#endif
}


UNAP::lduMatrix::lduMatrix
(
    const label&       nCells,
    const labelField&  lowerAddr,
    const labelField&  upperAddr,
    const scalarField& lower,
    const scalarField& diag,
    const scalarField& upper,
    const bool         reUse
)
:
    nCells_(nCells),
    lowerAddrPtr_(NULL),
    upperAddrPtr_(NULL),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    interfacesPtr_(NULL),
    losortPtr_(NULL),
    ownerStartPtr_(NULL),
    losortStartPtr_(NULL)
{
    if(reUse)
    {
        lowerAddrPtr_ = const_cast<labelField*> (&lowerAddr);
        upperAddrPtr_ = const_cast<labelField*> (&upperAddr);
        diagPtr_      = const_cast<scalarField*> (&diag);
        lowerPtr_     = const_cast<scalarField*> (&lower);
        upperPtr_     = const_cast<scalarField*> (&upper);
    }
    else
    {
        ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelField)
        ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelField)
        ALLOCATE_POINTER(diagPtr_,      diag,      scalarField)
        ALLOCATE_POINTER(upperPtr_,     upper,     scalarField)

        if((&lower) == (&upper))
        {
            lowerPtr_ = upperPtr_;
        }
        else
        {
            ALLOCATE_POINTER(lowerPtr_, lower, scalarField)
        }
    }

#ifdef SW_SLAVE
    mlbIter_ = NULL;
    rssIter_ = NULL;
    unatIter_ = NULL;
#endif
}


UNAP::lduMatrix::lduMatrix
(
    const label&       nCells,
    const labelField&  lowerAddr,
    const labelField&  upperAddr
)
:
    nCells_(nCells),
    lowerAddrPtr_(NULL),
    upperAddrPtr_(NULL),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    interfacesPtr_(NULL),
    losortPtr_(NULL),
    ownerStartPtr_(NULL),
    losortStartPtr_(NULL)
{
    ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelField)
    ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelField)

#ifdef SW_SLAVE
    mlbIter_ = NULL;
    rssIter_ = NULL;
    unatIter_ = NULL;
#endif
}


UNAP::lduMatrix::lduMatrix
(
    const label&       nCells,
    const labelField&  lowerAddr,
    const labelField&  upperAddr,
    const bool         reUse
)
:
    nCells_(nCells),
    lowerAddrPtr_(NULL),
    upperAddrPtr_(NULL),
    lowerPtr_(NULL),
    diagPtr_(NULL),
    upperPtr_(NULL),
    interfacesPtr_(NULL),
    losortPtr_(NULL),
    ownerStartPtr_(NULL),
    losortStartPtr_(NULL)
{
    if(reUse)
    {
        lowerAddrPtr_ = const_cast<labelField*> (&lowerAddr);
        upperAddrPtr_ = const_cast<labelField*> (&upperAddr);
    }
    else
    {
        ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelField)
        ALLOCATE_POINTER(upperAddrPtr_, upperAddr, labelField)
    }

#ifdef SW_SLAVE
    mlbIter_ = NULL;
    rssIter_ = NULL;
    unatIter_ = NULL;
#endif
}


UNAP::lduMatrix::~lduMatrix()
{
    if(lowerPtr_ != upperPtr_)
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

UNAP::labelField& UNAP::lduMatrix::lowerAddr() const
{
#ifdef DEBUG
	CHECK_POINTER(lowerAddrPtr_)
#endif

    return *lowerAddrPtr_;
}

UNAP::labelField& UNAP::lduMatrix::upperAddr() const
{
#ifdef DEBUG
    CHECK_POINTER(upperAddrPtr_)
#endif

    return *upperAddrPtr_;
}

UNAP::scalarField& UNAP::lduMatrix::lower() const
{
#ifdef DEBUG
    CHECK_POINTER(lowerPtr_)
#endif

    return *lowerPtr_;
}


UNAP::scalarField& UNAP::lduMatrix::diag() const
{
#ifdef DEBUG
    CHECK_POINTER(diagPtr_)
#endif

    return *diagPtr_;
}


UNAP::scalarField& UNAP::lduMatrix::upper() const
{
#ifdef DEBUG
    CHECK_POINTER(upperPtr_)
#endif

    return *upperPtr_;
}


void UNAP::lduMatrix::calcOwnerStart() const
{
    const labelField& own = lowerAddr();

    ownerStartPtr_ = new labelField(nCells() + 1, own.size());

    labelField& ownStart = *ownerStartPtr_;

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


const UNAP::labelField& UNAP::lduMatrix::ownerStartAddr() const
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
    labelField nNbrOfFace(nCells(), 0);

    const labelField& nbr = upperAddr();

    forAll(nbrI, nbr.size())
    {
        nNbrOfFace[nbr[nbrI]]++;
    }

    //- create temporary neighbor addressing
    label** cellNbrFaces = new label* [nCells_];
    forAll(cellI, nCells())
    {
        cellNbrFaces[cellI] = new label [nNbrOfFace[cellI]];
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
    losortPtr_ = new labelField(nbr.size(), -1);

    labelField& lst = *losortPtr_;

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
        delete [] cellNbrFaces[cellI];
        cellNbrFaces[cellI] = NULL;
    }
    delete [] cellNbrFaces;
    cellNbrFaces = NULL;
}


const UNAP::labelField& UNAP::lduMatrix::losortAddr() const
{
    if (!losortPtr_)
    {
        calcLosort();
    }

    return *losortPtr_;
}


void UNAP::lduMatrix::calcLosortStart() const
{
    losortStartPtr_ = new labelField(nCells() + 1, 0);

    labelField& lsrtStart = *losortStartPtr_;

    const labelField& nbr = upperAddr();

    const labelField& lsrt = losortAddr();

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


const UNAP::labelField& UNAP::lduMatrix::losortStartAddr() const
{
    if (!losortStartPtr_)
    {
        calcLosortStart();
    }

    return *losortStartPtr_;
}


void UNAP::lduMatrix::createInterfacesTopology
(
    const label   nNeiProcs,
    const label*  destRank,
    const label*  offDiagRows,
    const label*  offDiagStarts
)
{
    PtrList<patch>* patchesPtr = new PtrList<patch>(nNeiProcs);

    forAll(intI, nNeiProcs)
    {
        const label neighborID = destRank[intI] - 1;
        const label localSize = offDiagStarts[intI+1] - offDiagStarts[intI];
        patch* patchIPtr = new patch(localSize, MYID, neighborID);

        label* localFaceCells = new label[localSize];

        forAll(faceI, localSize)
        {
            label start = offDiagStarts[intI] + faceI - 1;
            localFaceCells[faceI] = offDiagRows[start] - 1;
        }

        labelField* locFaceCellsPtr = new labelField(localFaceCells, localSize, true);
        scalarField* localDataPtr = new scalarField(localSize);

        patchIPtr->faceCells(*locFaceCellsPtr);
        patchIPtr->patchCoeffs(*localDataPtr);

        patchesPtr->setLevel(intI, *patchIPtr);
    }

    interfaces* interfacesLocalPtr = new interfaces(*patchesPtr);
    matrixInterfaces(*interfacesLocalPtr);
}

void UNAP::lduMatrix::fillInterfacesCofficients
(
    const label*   offDiagStarts,
    const scalar*  offDiagCoeffs
)
{
    interfaces& lduAInter = matrixInterfaces();

    const label nNeiProcs = lduAInter.size();

    forAll(intI, nNeiProcs)
    {
        patch& patchI = lduAInter.patchList(intI);
        const label localSize = offDiagStarts[intI+1] - offDiagStarts[intI];
        scalar* localData = patchI.patchCoeffs().begin();

        forAll(faceI, localSize)
        {
            label start = offDiagStarts[intI] + faceI - 1;
            localData[faceI] = offDiagCoeffs[start];
        }
    }
}


void UNAP::lduMatrix::initInterfaces(const scalarField& psi) const
{
    if(PARRUN)
    {
        matrixInterfaces().initMatrixInterfaces(psi);
    }
}


void UNAP::lduMatrix::updateInterfaces(scalarField& Apsi) const
{
    if(PARRUN)
    {
        matrixInterfaces().updateMatrixInterfaces(Apsi);
    }
}


void UNAP::lduMatrix::setMatrixCoeffients
(
    const scalarField& diag,
    const scalarField& upper,
    const scalarField& lower
)
{
    nCells_ = diag.size();
    ALLOCATE_POINTER(diagPtr_,  diag,  scalarField)

    ALLOCATE_POINTER(upperPtr_, upper, scalarField)

    if((&lower) == (&upper))
    {
        lowerPtr_ = upperPtr_;
    }
    else
    {
        ALLOCATE_POINTER(lowerPtr_, lower, scalarField)
    }
}


void UNAP::lduMatrix::setMatrixCoeffients
(
    const scalarField& diag,
    const scalarField& upper,
    const scalarField& lower,
    const bool reUse
)
{
    nCells_ = diag.size();
    if(reUse)
    {
        DELETE_OBJECT_POINTER(diagPtr_)
        this->diagPtr_  = const_cast<scalarField*>(&diag);

        DELETE_OBJECT_POINTER(upperPtr_)
        this->upperPtr_ = const_cast<scalarField*>(&upper);

        if((&lower) == (&upper))
        {
            this->lowerPtr_ = this->upperPtr_;
        }
        else
        {
            DELETE_OBJECT_POINTER(lowerPtr_)
            this->lowerPtr_ = const_cast<scalarField*>(&lower);
        }
    }
    else
    {
        setMatrixCoeffients(diag, upper, lower);
    }
}



void UNAP::lduMatrix::setMatrixTopology
(
    const labelField& upperAddr,
    const labelField& lowerAddr
)
{
    ALLOCATE_POINTER(upperAddrPtr_,  upperAddr,  labelField)

    ALLOCATE_POINTER(lowerAddrPtr_, lowerAddr, labelField)
}


void UNAP::lduMatrix::setMatrixTopology
(
    const labelField& upperAddr,
    const labelField& lowerAddr,
    const bool reUse
)
{
    if(reUse)
    {
        DELETE_OBJECT_POINTER(upperAddrPtr_)
        this->upperAddrPtr_  = const_cast<labelField*>(&upperAddr);

        DELETE_OBJECT_POINTER(lowerAddrPtr_)
        this->lowerAddrPtr_  = const_cast<labelField*>(&lowerAddr);
    }
    else
    {
        setMatrixTopology(upperAddr, lowerAddr);
    }
}

#ifdef SW_SLAVE

void UNAP::lduMatrix::constructMLBIterator()
{
    label* uPtr = this->upperAddr().begin();
    label* lPtr = this->lowerAddr().begin();
    label nFaces = this->lowerAddr().size();
    label nCells = this->size();

    if(nFaces >= SpMVAccSize)
    {
        label* cellWeights = new label[nCells];
        label* faceWeights = new label[nFaces];

        for(int i=0;i<nCells;i++) {cellWeights[i] = 3;}
        for(int i=0;i<nFaces;i++) {faceWeights[i] = 2;}
        mlbIter_ = constructMLBIteratorFromUNAP(lPtr, uPtr,
                    cellWeights, faceWeights, nFaces);

        this->unatEdgeMap_ = getEdgeMap(this->mlbIter_);
        this->unatCellMap_ = getCellMap(this->mlbIter_);

        this->mlbIter_->reorderEdges(lPtr, uPtr, nFaces, nCells);
        if(PARRUN)
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
    scalar* uPtr = this->upper().begin();
    scalar* lPtr = this->lower().begin();
    scalar* dPtr = this->diag().begin();
    label nFaces = this->lower().size();
    label nCells = this->size();

    if(mlbIter_)
    {
        Arrays selfConnData, frontEdgeData, backEdgeData;
        constructSingleArray(selfConnData, 1, nCells, COPYIN,
                        (swFloat*)dPtr);
        constructSingleArray(frontEdgeData, 1, nFaces, COPYIN,
                        (swFloat*)uPtr);
        constructSingleArray(backEdgeData, 1, nFaces, COPYIN,
                        (swFloat*)lPtr);
        this->mlbIter_->reorderVertexData(&selfConnData);
        this->mlbIter_->reorderEdgeData(&backEdgeData, &frontEdgeData);
    }
}


void UNAP::lduMatrix::restoreVector(scalarField& vv)
{
    scalar* vvPtr = vv.begin();
    label nCells = vv.size();

    if(mlbIter_)
    {
        Arrays selfConnData;
        constructSingleArray(selfConnData, 1, nCells, COPYIN,
                        (swFloat*)vvPtr);
        this->mlbIter_->restoreVertexData(&selfConnData);
    }
}


void UNAP::lduMatrix::reorderVector(scalarField& vv)
{
    scalar* vvPtr = vv.begin();
    label nCells = vv.size();

    if(mlbIter_)
    {
        Arrays selfConnData;
        constructSingleArray(selfConnData, 1, nCells, COPYIN,
                        (swFloat*)vvPtr);
        this->mlbIter_->reorderVertexData(&selfConnData);
    }
}



void UNAP::lduMatrix::constructRSSIterator()
{
    label* uPtr = this->upperAddr().begin();
    label* lPtr = this->lowerAddr().begin();
    label nFaces = this->lowerAddr().size();
    label nCells = this->size();

    if(nFaces >= SpMVAccSize)
    {
        label* cellWeights = new label[nCells];
        label* faceWeights = new label[nFaces];

        for(int i=0;i<nCells;i++) {cellWeights[i] = 3;}
        for(int i=0;i<nFaces;i++) {faceWeights[i] = 2;}
        rssIter_ = constructRSSIteratorFromUNAP(lPtr, uPtr,
                    cellWeights, faceWeights, nFaces);

        DELETE_POINTER(cellWeights)
        DELETE_POINTER(faceWeights)
    }

    unatIter_ = rssIter_;
}


#endif

