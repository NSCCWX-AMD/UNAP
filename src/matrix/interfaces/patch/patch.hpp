#ifndef PATCH_HPP
#define PATCH_HPP

#include "unapVector.hpp"

namespace UNAP
{
class patch
{
private:
  //- number of faces in this patch
  label size_;

  //- current processor number
  label myProcNo_;

  //- neighbor processor number
  label neighbProcNo_;

  //- face-cell addressing
  mutable labelVector *faceCellsPtr_;

  mutable labelVector *faceCells2Ptr_;

  //- coefficients
  mutable scalarVector *patchCoeffsPtr_;

  //- face restriction addressing if needed
  mutable labelVector *faceRestrictAddressingPtr_;

public:
  patch(label size, label myProcNo, label neighbProcNo);

  virtual ~patch();

  //- return neighbor processor number
  inline label neighbProcNo() const { return neighbProcNo_; }

  inline void neighbProcNo(const label i) { neighbProcNo_ = i; }

  //- return current processor number
  inline label myProcNo() const { return myProcNo_; }

  inline void myProcNo(const label i) { myProcNo_ = i; }

  //- return faceCells
  inline labelVector &faceCells() const
  {
    CHECK_POINTER(faceCellsPtr_)
    return *faceCellsPtr_;
  }

  inline labelVector &faceCells2() const
  {
    CHECK_POINTER(faceCells2Ptr_)
    return *faceCells2Ptr_;
  }

  inline void faceCells(labelVector &a) const
  {
    DELETE_OBJECT_POINTER(faceCellsPtr_)
    faceCellsPtr_ = &a;
  }

  inline void faceCells2(labelVector &a) const
  {
    DELETE_OBJECT_POINTER(faceCells2Ptr_)
    faceCells2Ptr_ = &a;
  }

  //- return patch coefficients
  inline scalar patchCoeffs(const label faceI) const
  {
    return (*patchCoeffsPtr_)[faceI];
  }

  inline scalarVector &patchCoeffs() const
  {
    CHECK_POINTER(patchCoeffsPtr_)
    return *patchCoeffsPtr_;
  }

  inline void patchCoeffs(scalarVector &a)
  {
    DELETE_OBJECT_POINTER(patchCoeffsPtr_)
    patchCoeffsPtr_ = &a;
  }

  inline label size() const { return size_; }

  inline void size(const label i) { size_ = i; }

  inline void faceRestrictAddressing(labelVector &a) const
  {
    DELETE_OBJECT_POINTER(faceRestrictAddressingPtr_)
    faceRestrictAddressingPtr_ = &a;
  }

  inline labelVector &faceRestrictAddressing() const
  {
    CHECK_POINTER(faceRestrictAddressingPtr_)
    return *faceRestrictAddressingPtr_;
  }

  void reorderPatchFaceCells(const label *cellMap);
};

}  // namespace UNAP
#endif  //- PATCH_HPP
