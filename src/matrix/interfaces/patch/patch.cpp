#include "patch.hpp"

UNAP::patch::patch(label size, label myProcNo, label neighbProcNo)
    : size_(size),
      myProcNo_(myProcNo),
      neighbProcNo_(neighbProcNo),
      faceCellsPtr_(NULL),
      patchCoeffsPtr_(NULL),
      faceRestrictAddressingPtr_(NULL)
{
}

UNAP::patch::~patch()
{
  DELETE_OBJECT_POINTER(faceCellsPtr_)
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
}
