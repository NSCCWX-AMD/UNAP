#include "matrixConversion.hpp"

UNAP::lduMatrix &UNAP::coo2ldu(const scalar *dataPtr,
                               const label *rowsPtr,
                               const label *columnPtr,
                               const label nCells,
                               const label size,
                               const bool symm,
                               Communicator *other_comm)
{
  //- input coo matrix must have been arranged in order of row
  //- in each row, data must have been arranged in order of column
  //- structure of matrix must be symmetric
  //- otherwise this part will not work and the results should be wrong
  if (!other_comm)
  {
    std::cout << "Error: " << __FILE__ << " in " << __LINE__
              << "The Communicator is NULL !" << ENDL;
    ERROR_EXIT;
  }
  //- number of no-zero in upper
  const label nZeros = (size - nCells) / 2;

  scalarVector diag(nCells, other_comm);
  scalarVector upper(nZeros, other_comm);

  labelVector upperAddr(nZeros, other_comm);
  labelVector lowerAddr(nZeros, other_comm);

  //- asymmetric case
  scalarVector *lowerPtr = NULL;
  labelVector *lowerColPtr = NULL;
  labelVector *upperNbrsPtr = NULL;

  if (!symm)
  {
    lowerPtr = new scalarVector(nZeros, other_comm);
    lowerColPtr = new labelVector(nZeros, other_comm);
    upperNbrsPtr = new labelVector(nCells, 0, other_comm);
  }
  else
  {
    lowerPtr = &upper;
  }

  scalarVector &lower = *lowerPtr;
  labelVector &lowerCol = *lowerColPtr;
  labelVector &upperNbrs = *upperNbrsPtr;

  label upperCount = 0;
  label lowerCount = 0;

  forAll(i, size)
  {
    label row = rowsPtr[i];
    label col = columnPtr[i];

    if (row == col)
    {
      diag[row] = dataPtr[i];
    }
    else if (col > row)
    {
      upperAddr[upperCount] = col;
      lowerAddr[upperCount] = row;
      upper[upperCount] = dataPtr[i];
      upperCount++;

      if (!symm)
      {
        upperNbrs[row]++;
      }
    }
    else if ((col < row) && (!symm))
    {
      lowerCol[lowerCount] = col;
      lower[lowerCount] = dataPtr[i];
      lowerCount++;
    }
  }

  //- sort lower data
  if (!symm)
  {
    sortData(lower, lowerCol, upperNbrs, other_comm);
  }

#ifdef DEBUG
  if (upperCount != nZeros)
  {
    other_comm->log()
        << "ERROR in " << __FILE__ << " " << __LINE__
        << ": the input COO matrix is not a structural symmetric matrix!"
        << ENDL;
    ERROR_EXIT;
  }
#endif

  lduMatrix *lduAPtr = new lduMatrix(
      nCells, lowerAddr, upperAddr, lower, diag, upper, other_comm);

  if (!symm)
  {
    DELETE_OBJECT_POINTER(lowerPtr)
    DELETE_OBJECT_POINTER(lowerColPtr)
    DELETE_OBJECT_POINTER(upperNbrsPtr)
  }

  return *lduAPtr;
}

UNAP::lduMatrix &UNAP::csr2ldu(const scalar *dataPtr,
                               const label *compRowsPtr,
                               const label *columnPtr,
                               const label nCells,
                               const label size,
                               const bool symm,
                               Communicator *other_comm)
{
  label *rowsPtr = new label[size];

  forAll(i, nCells)
  {
    label nnzInRow = compRowsPtr[i + 1] - compRowsPtr[i];
    forAll(j, nnzInRow) { rowsPtr[compRowsPtr[i] + j] = i; }
  }

  lduMatrix &lduA =
      coo2ldu(dataPtr, rowsPtr, columnPtr, nCells, size, symm, other_comm);
  delete[] rowsPtr;
  return lduA;
}
