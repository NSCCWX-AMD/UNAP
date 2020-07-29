#ifndef LDUAGGLOMERATION_HPP
#define LDUAGGLOMERATION_HPP

#include "lduMatrix.hpp"

namespace UNAP
{
class lduAgglomeration : public matrix::agglomeration
{
private:
  //- direction of cell loop for the current level
  static bool forward_;

  //- max number of levels
  label maxLevels_;

  //- number of cells in coarsest level
  label nCellsInCoarsestLevel_;

  //- the number of cells in each level
  labelVector nCells_;

  //- number of levels created
  label nCreatedLevels_;

  //- number of levels to merge, 1 = don't merge, 2 = merge pairs etc.
  label mergeLevels_;

  //- finest matrix
  const lduMatrix &finestMatrix_;

  //- hierarchy of coarse matrix levels
  PtrList<lduMatrix> coarseMatrixLevels_;

  //- cell restriction addressing array.
  //- maps from the finer to the coarser level.
  PtrList<labelVector> restrictAddressing_;

  //- face restriction addressing array.
  //- maps from the finer to the coarser level.
  //- positive indices map the finer faces which form part of the boundary
  //  of the coarser cells to the corresponding coarser cell face.
  //- negative indices map the finer faces which are internal to the
  //  coarser cells to minus the corresponding coarser cell index minus 1.
  PtrList<labelVector> faceRestrictAddressing_;

  //- face restriction addressing in interface patch
  //- maps from the finer to the coarser level.
  PtrList<labelVector> patchFaceRestrictAddressing_;

  //- assemble coarse mesh addressing
  void agglomerateLduAddressing(const label fineLevelIndex);

  //- shrink the number of levels to that specified
  void compactLevels(const label nCreatedLevels);

  //- check the need for further agglomeration
  bool continueAgglomerating(const label nCoarseCells) const;

  //- combine levels
  void combineLevels(const label curLevel);

  //- calculate and return agglomeration of given level
  labelVector &agglomerate(label &nCoarseCells,
                           const lduMatrix &fineA,
                           const scalarVector &weights);

  //- agglomerate coarse matrix
  virtual void agglomerateMatrix(const label fineLevelIndex);

  const lduMatrix &matrixLevel(const label leveli) const;

public:
  //- constructors
  lduAgglomeration(const lduMatrix &A);

  //- destructor
  virtual ~lduAgglomeration();

  //- agglomerate all levels starting from the given face weights
  virtual void agglomerate(const scalarVector &weights);

  virtual label size() const { return nCreatedLevels_; }

  //- access to coarse matrix level
  virtual const lduMatrix &coarseMatrixLevels(const label leveli) const
  {
    return coarseMatrixLevels_[leveli];
  }

  lduMatrix &coarseMatrix(const label leveli) const
  {
    return coarseMatrixLevels_[leveli];
  }

  //- access to restrictAddressing_
  //- return cell restrict addressing of given level
  virtual const labelVector &restrictAddressing(const label leveli) const
  {
    return restrictAddressing_[leveli];
  }

  //- access to faceRestrictAddressing_
  //- return face restrict addressing of given level
  virtual const labelVector &faceRestrictAddressing(const label leveli) const
  {
    return faceRestrictAddressing_[leveli];
  }

  virtual void SET_nCellsInCoarsestLevel(const label i);

  virtual void SET_maxLevels(const label i);

  void agglomerationReorderTopo();

  void agglConstructRSSIterator();
};

}  // namespace UNAP

#endif  //- LDUAGGLOMERATION_HPP
