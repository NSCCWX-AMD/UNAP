#include "patch.hpp"

UNAP::patch::patch(label size, label myProcNo, label neighbProcNo)
    : size_(size),
      myProcNo_(myProcNo),
      neighbProcNo_(neighbProcNo),
      faceCellsPtr_(NULL),
      faceCells2Ptr_(NULL),
      patchCoeffsPtr_(NULL),
      faceRestrictAddressingPtr_(NULL)
{
}

UNAP::patch::~patch()
{
  DELETE_OBJECT_POINTER(faceCellsPtr_)
  DELETE_OBJECT_POINTER(faceCells2Ptr_)
  DELETE_OBJECT_POINTER(patchCoeffsPtr_)
  DELETE_OBJECT_POINTER(faceRestrictAddressingPtr_)
}

void UNAP::patch::reorderPatchFaceCells(const label *cellMap)
{
  forAll(faceI, size_)
  {
    label *faceCellsPtr = faceCellsPtr_->begin();
    faceCellsPtr[faceI] = cellMap[faceCellsPtr[faceI]];
  }

  if (myProcNo_ == neighbProcNo_)
  {
    forAll(faceI, size_)
    {
      label *faceCells2Ptr = faceCells2Ptr_->begin();
      faceCells2Ptr[faceI] = cellMap[faceCells2Ptr[faceI]];
    }
  }
}
