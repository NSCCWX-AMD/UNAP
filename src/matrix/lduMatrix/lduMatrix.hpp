#ifndef LDUMATRIX_HPP
#define LDUMATRIX_HPP

#include "unapMatrix.hpp"

#ifdef SW_SLAVE
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"
#include "iterator.hpp"
#endif

//- to do
//- add API for reusing data

namespace UNAP
{
class lduMatrix
:
	public matrix
{
private:

    //- number of cells
    label nCells_;

    //- addressing
    labelField *lowerAddrPtr_, *upperAddrPtr_;

	//- coefficients (not including interfaces)
    scalarField *lowerPtr_, *diagPtr_, *upperPtr_;

    //- interfaces
    interfaces* interfacesPtr_;

    //- losort addressing
    mutable labelField* losortPtr_;

    //- owner start addressing
    mutable labelField* ownerStartPtr_;

    //- losort start addressing
    mutable labelField* losortStartPtr_;

    //- calculate losort
    void calcLosort() const;

    //- calculate owner start
    void calcOwnerStart() const;

    //- calculate losort start
    void calcLosortStart() const;

#ifdef SW_SLAVE
    //- multi-level block iterator
    UNAT::MultiLevelBlockIterator *mlbIter_;

    //- start from 1
    //- has negative value, which stands for lower
    label* unatEdgeMap_;

    //- start from 0
    label* unatCellMap_;


    //- row subsection iterator
    UNAT::RowSubsectionIterator* rssIter_;


    UNAT::Iterator* unatIter_;


#endif


public:
    //- constructors

    lduMatrix();

    lduMatrix
    (
        const label& nCells,
        const labelField& lowerAddr,
        const labelField& upperAddr,
        const scalarField& lower,
        const scalarField& diag,
        const scalarField& upper
    );

    lduMatrix
    (
        const label& nCells,
        const labelField& lowerAddr,
        const labelField& upperAddr,
        const scalarField& lower,
        const scalarField& diag,
        const scalarField& upper,
        const bool reUse
    );

    //- constructor
    //- only topology
    lduMatrix
    (
        const label& nCells,
        const labelField& lowerAddr,
        const labelField& upperAddr
    );

    lduMatrix
    (
        const label& nCells,
        const labelField& lowerAddr,
        const labelField& upperAddr,
        const bool reUse
    );

    //- destructor
    ~lduMatrix();


    //- access to addressing
    virtual labelField& lowerAddr() const;
    virtual labelField& upperAddr() const;

	//- access to coefficients
    virtual scalarField& lower() const;
    virtual scalarField& diag () const;
    virtual scalarField& upper() const;

    void SET_lowerAddr(labelField& newLowerAddr)
    {
        ALLOCATE_POINTER(lowerAddrPtr_, newLowerAddr, labelField)
    }

    void SET_upperAddr(labelField& newUpperAddr)
    {
        ALLOCATE_POINTER(upperAddrPtr_, newUpperAddr, labelField)
    }

    void SET_lower(scalarField& newLower)
    {
        if(this->symm())
        {
            lowerPtr_ = NULL;
            // COUT << "Here works" << ENDL;
        }

        ALLOCATE_POINTER(lowerPtr_, newLower, scalarField)
    }

    void SET_lower(label newSize)
    {
        DELETE_OBJECT_POINTER(lowerPtr_)

        lowerPtr_ =  new scalarField(newSize);
    }

    void SET_upper(scalarField& newUpper)
    {
        ALLOCATE_POINTER(upperPtr_, newUpper, scalarField)
    }

    void SET_upper(label newSize)
    {
        DELETE_OBJECT_POINTER(upperPtr_)

        upperPtr_ =  new scalarField(newSize);
    }

    void SET_diag(scalarField& newDiag)
    {
        ALLOCATE_POINTER(diagPtr_, newDiag, scalarField)
    }

    void SET_diag(label newSize)
    {
        DELETE_OBJECT_POINTER(diagPtr_)

        diagPtr_ =  new scalarField(newSize);
    }

    scalarField& diag()
    {
        return *diagPtr_;
    }

    void freeLower()
    {
        if(!this->symm())
        {
            DELETE_OBJECT_POINTER(lowerPtr_);
            lowerPtr_ = upperPtr_;
        }
    }

    void setSymm()
    {
        if(upperPtr_)
        {
            lowerPtr_ = upperPtr_;
        }
        else if(lowerPtr_)
        {
            upperPtr_ = lowerPtr_;
        }
    }

    //- access to symmetric or asymmetric
    virtual bool symm() const
    {
        if(lowerPtr_ == upperPtr_)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    virtual label size() const
    {
        return nCells_;
    }

    //- access to nCells
    virtual label nCells() const
    {
        return nCells_;
    }

    void size(label size)
    {
        nCells_ = size;
    }


    virtual void spMV
    (
    	scalarField& Apsi,
		const scalarField& psi
    ) const;

    //- return losort addressing
    const labelField& losortAddr() const;

    //- return owner start addressing
    const labelField& ownerStartAddr() const;

    //- return losort start addressing
    const labelField& losortStartAddr() const;

    virtual interfaces& matrixInterfaces() const
    {
        return *interfacesPtr_;
    }

    virtual void matrixInterfaces(interfaces& a)
    {
        interfacesPtr_ = &a;
    }

    void matrixInterfaces(const label size)
    {
        PtrList<patch>* patchesPtr = new PtrList<patch> (size);
        interfacesPtr_ = new interfaces(*patchesPtr);
    }

    //- initialize interfaces
    virtual void initInterfaces(const scalarField& psi) const;

    //- update interfaces
    virtual void updateInterfaces(scalarField& Apsi) const;

    //- fill coefficients
    void setMatrixCoeffients
    (
        const scalarField& diag,
        const scalarField& upper,
        const scalarField& lower
    );

    void setMatrixCoeffients
    (
        const scalarField& diag,
        const scalarField& upper,
        const scalarField& lower,
        const bool reuse
    );

    //- fill topology
    void setMatrixTopology
    (
        const labelField& upperAddr,
        const labelField& lowerAddr,
        const bool reUse
    );

    void setMatrixTopology
    (
        const labelField& upperAddr,
        const labelField& lowerAddr
    );

#ifdef SW_SLAVE
    void constructMLBIterator();

    UNAT::MultiLevelBlockIterator* mlbIter() const
    {
        return mlbIter_;
    }

    void reorderLDUValues();

    void restoreVector(scalarField& vv);

    void reorderVector(scalarField& vv);

    const label* unatEdgeMap() const
    {
        return unatEdgeMap_;
    }

    const label* unatCellMap() const
    {
        return unatCellMap_;
    }


    void constructRSSIterator();


#endif
};

} //- namespace UNAP
#endif //- LDUMATRIX_HPP
