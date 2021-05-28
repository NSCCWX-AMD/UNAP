#ifndef LDUMATRIX_HPP
#define LDUMATRIX_HPP

#include "unapMatrix.hpp"

#ifdef SW_SLAVE
#include "iterator.hpp"
#include "multiLevelBlockIterator.hpp"
#include "rowSubsectionIterator.hpp"
#endif

//- to do
//- add API for reusing data

namespace UNAP
{
class lduMatrix : public matrix
{
 private:
  //- number of cells
  label nCells_;

  //- addressing
  labelVector *lowerAddrPtr_, *upperAddrPtr_;

  //- coefficients (not including interfaces)
  scalarVector *lowerPtr_, *diagPtr_, *upperPtr_;

  //- interfaces
  interfaces *interfacesPtr_;

  //- losort addressing
  mutable labelVector *losortPtr_;

  //- owner start addressing
  mutable labelVector *ownerStartPtr_;

  //- losort start addressing
  mutable labelVector *losortStartPtr_;

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
  label *unatEdgeMap_;

  //- start from 0
  label *unatCellMap_;

  //- row subsection iterator
  UNAT::RowSubsectionIterator *rssIter_;

  UNAT::Iterator *unatIter_;

#endif

 public:
  //- constructors

  lduMatrix(Communicator *other_comm);

  lduMatrix(const label &nCells,
            const labelVector &lowerAddr,
            const labelVector &upperAddr,
            const scalarVector &lower,
            const scalarVector &diag,
            const scalarVector &upper,
            Communicator *other_comm);

  lduMatrix(const label &nCells,
            const labelVector &lowerAddr,
            const labelVector &upperAddr,
            const scalarVector &lower,
            const scalarVector &diag,
            const scalarVector &upper,
            const bool reUse,
            Communicator *other_comm);

  //- constructor
  //- only topology
  lduMatrix(const label &nCells,
            const labelVector &lowerAddr,
            const labelVector &upperAddr,
            Communicator *other_comm);

  lduMatrix(const label &nCells,
            const labelVector &lowerAddr,
            const labelVector &upperAddr,
            const bool reUse,
            Communicator *other_comm);

  //- destructor
  ~lduMatrix();

  //- access to addressing
  virtual labelVector &lowerAddr() const;
  virtual labelVector &upperAddr() const;

  //- access to coefficients
  virtual scalarVector &lower() const;
  virtual scalarVector &diag() const;
  virtual scalarVector &upper() const;

  void SET_lowerAddr(labelVector &newLowerAddr)
  {
    ALLOCATE_POINTER(lowerAddrPtr_, newLowerAddr, labelVector)
  }

  void SET_upperAddr(labelVector &newUpperAddr)
  {
    ALLOCATE_POINTER(upperAddrPtr_, newUpperAddr, labelVector)
  }

  void SET_lower(scalarVector &newLower)
  {
    if (this->symm())
    {
      lowerPtr_ = NULL;
      // COUT << "Here works" << ENDL;
    }

    ALLOCATE_POINTER(lowerPtr_, newLower, scalarVector)
  }

  void SET_lower(label newSize)
  {
    DELETE_OBJECT_POINTER(lowerPtr_)

    lowerPtr_ = new scalarVector(newSize, this->commcator_);
  }

  void SET_upper(scalarVector &newUpper)
  {
    ALLOCATE_POINTER(upperPtr_, newUpper, scalarVector)
  }

  void SET_upper(label newSize)
  {
    DELETE_OBJECT_POINTER(upperPtr_)

    upperPtr_ = new scalarVector(newSize, this->commcator_);
  }

  void SET_diag(scalarVector &newDiag)
  {
    ALLOCATE_POINTER(diagPtr_, newDiag, scalarVector)
  }

  void SET_diag(label newSize)
  {
    DELETE_OBJECT_POINTER(diagPtr_)

    diagPtr_ = new scalarVector(newSize, this->commcator_);
  }

  scalarVector &diag() { return *diagPtr_; }

  void freeLower()
  {
    if (!this->symm())
    {
      DELETE_OBJECT_POINTER(lowerPtr_);
      lowerPtr_ = upperPtr_;
    }
  }

  void setSymm()
  {
    if (upperPtr_)
    {
      lowerPtr_ = upperPtr_;
    }
    else if (lowerPtr_)
    {
      upperPtr_ = lowerPtr_;
    }
  }

  //- access to symmetric or asymmetric
  virtual bool symm() const
  {
    if (lowerPtr_ == upperPtr_)
    {
      return true;
    }
    else
    {
      return false;
    }
  }

  virtual label size() const { return nCells_; }

  //- access to nCells
  virtual label nCells() const { return nCells_; }

  void size(label size) { nCells_ = size; }

  virtual void spMV(scalarVector &Apsi, const scalarVector &psi) const;

  //- return losort addressing
  const labelVector &losortAddr() const;

  //- return owner start addressing
  const labelVector &ownerStartAddr() const;

  //- return losort start addressing
  const labelVector &losortStartAddr() const;

  virtual interfaces &matrixInterfaces() const { return *interfacesPtr_; }

  virtual void matrixInterfaces(interfaces &a) { interfacesPtr_ = &a; }

  void matrixInterfaces(const label size)
  {
    PtrList<patch> *patchesPtr = new PtrList<patch>(size);
    interfacesPtr_ = new interfaces(*patchesPtr, this->commcator_);
  }

  void createInterfacesTopology(const label nNeiProcs,
                                const label *destRank,
                                const label *offDiagRows,
                                const label *offDiagStarts);

  void fillInterfacesCofficients(const label *offDiagStarts,
                                 const scalar *offDiagCoeffs);

  //- initialize interfaces
  virtual void initInterfaces(const scalarVector &psi) const;

  //- update interfaces
  virtual void updateInterfaces(scalarVector &Apsi) const;

  //- fill coefficients
  void setMatrixCoeffients(const scalarVector &diag,
                           const scalarVector &upper,
                           const scalarVector &lower);

  void setMatrixCoeffients(const scalarVector &diag,
                           const scalarVector &upper,
                           const scalarVector &lower,
                           const bool reuse);

  //- fill topology
  void setMatrixTopology(const labelVector &upperAddr,
                         const labelVector &lowerAddr,
                         const bool reUse);

  void setMatrixTopology(const labelVector &upperAddr,
                         const labelVector &lowerAddr);

#ifdef SW_SLAVE
  void constructMLBIterator();

  UNAT::MultiLevelBlockIterator *mlbIter() const { return mlbIter_; }

  void reorderLDUValues();

  void restoreVector(scalarVector &vv);

  void reorderVector(scalarVector &vv);

  const label *unatEdgeMap() const { return unatEdgeMap_; }

  const label *unatCellMap() const { return unatCellMap_; }

  void setMlbIter(UNAT::MultiLevelBlockIterator *mlbIter);
  // TODO:temperary
  void unatMapFree()
  {
    DELETE_POINTER(unatEdgeMap_);
    DELETE_POINTER(unatCellMap_);
  }

  void constructRSSIterator();

#endif
};

}  // namespace UNAP
#endif  //- LDUMATRIX_HPP
