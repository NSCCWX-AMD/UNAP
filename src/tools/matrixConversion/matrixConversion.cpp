#include "matrixConversion.hpp"

#include "printUNAP.hpp"

void UNAP::sortData(scalarVector &data,
                    const labelVector &order,
                    const labelVector &cellFaces,
                    Communicator *other_comm)
{
  if (data.getCommunicator() != order.getCommunicator() &&
      cellFaces.getCommunicator() != data.getCommunicator())
  {
    data.getCommunicator()->log()
        << "Error: " << __FILE__ << " in " << __LINE__
        << "The Communicator between data , order and cellFaces are "
           "different!\n";
  }
  const label nCells = cellFaces.size();

  labelVector cellFaceOffsets(nCells + 1, other_comm);

  cellFaceOffsets[0] = 0;
  forAll(celli, nCells)
  {
    cellFaceOffsets[celli + 1] = cellFaceOffsets[celli] + cellFaces[celli];
  }

  //- size of faces in upper
  const label nFaces = data.size();

  scalarVector dataOld(data);

  //- count the data input
  labelVector count(nCells, 0, other_comm);

  forAll(facei, nFaces)
  {
    label cellLoc = order[facei];
    data[cellFaceOffsets[cellLoc] + count[cellLoc]] = dataOld[facei];
    count[cellLoc]++;
  }
}

void UNAP::reorderCOO(scalar *dataPtr,
                      label *rowsPtr,
                      label *colsPtr,
                      const label nCells,
                      const label size,
                      Communicator *other_comm)
{
  labelVector countsInRow(nCells, other_comm);

  forAll(i, size) { countsInRow[rowsPtr[i]]++; }

  labelVector countsInRowOffsets(nCells + 1, other_comm);

  forAll(i, nCells)
  {
    countsInRowOffsets[i + 1] = countsInRowOffsets[i] + countsInRow[i];
  }

  labelVector posInRow(size, other_comm);
  forAll(i, size)
  {
    label rowLocal = rowsPtr[i];
    label rowStartPos = countsInRowOffsets[rowLocal];
    label localRowSize =
        countsInRowOffsets[rowLocal + 1] - countsInRowOffsets[rowLocal];

    forAll(k, localRowSize)
    {
      if (colsPtr[i] > colsPtr[rowStartPos + k])
      {
        posInRow[i]++;
      }
    }
  }

  labelVector rowTemp(size, other_comm);
  labelVector colTemp(size, other_comm);
  scalarVector valTemp(size, other_comm);

  forAll(i, size)
  {
    label rowLocal = rowsPtr[i];
    label rowStartPos = countsInRowOffsets[rowLocal];
    label pos = rowStartPos + posInRow[i];

    valTemp[pos] = dataPtr[i];
    rowTemp[pos] = rowsPtr[i];
    colTemp[pos] = colsPtr[i];
  }

  forAll(i, size)
  {
    dataPtr[i] = valTemp[i];
    rowsPtr[i] = rowTemp[i];
    colsPtr[i] = colTemp[i];
  }
}

void UNAP::reorderValue(scalar *val,
                        const label *newOrder,
                        const label size,
                        Communicator *other_comm)
{
  scalarVector valTemp(val, size, false, other_comm);

  forAll(i, size) { val[newOrder[i]] = valTemp[i]; }
}

void UNAP::reorderUFace(label *rowsPtr,
                        label *colsPtr,
                        const label nCells,
                        const label size,
                        label *newOrder,
                        Communicator *other_comm)
{
  labelVector countsInRow(nCells, other_comm);

  forAll(i, size) { countsInRow[rowsPtr[i]]++; }

  labelVector countsInRowOffsets(nCells + 1, other_comm);

  forAll(i, nCells)
  {
    countsInRowOffsets[i + 1] = countsInRowOffsets[i] + countsInRow[i];
  }

  labelVector posInRow(size, other_comm);
  forAll(i, size)
  {
    label rowLocal = rowsPtr[i];
    label rowStartPos = countsInRowOffsets[rowLocal];
    label localRowSize =
        countsInRowOffsets[rowLocal + 1] - countsInRowOffsets[rowLocal];

    forAll(k, localRowSize)
    {
      if (colsPtr[i] > colsPtr[rowStartPos + k])
      {
        posInRow[i]++;
      }
    }
  }

  labelVector rowTemp(size, other_comm);
  labelVector colTemp(size, other_comm);

  forAll(i, size)
  {
    label rowLocal = rowsPtr[i];
    label rowStartPos = countsInRowOffsets[rowLocal];
    label pos = rowStartPos + posInRow[i];

    newOrder[i] = pos;
    rowTemp[pos] = rowsPtr[i];
    colTemp[pos] = colsPtr[i];
  }

  forAll(i, size)
  {
    rowsPtr[i] = rowTemp[i];
    colsPtr[i] = colTemp[i];
  }
}

void UNAP::reorderLFace(label *rowsPtr,
                        label *colsPtr,
                        const label nCells,
                        const label size,
                        label *newOrder,
                        Communicator *other_comm)
{
  labelVector countsInCol(nCells, other_comm);

  forAll(i, size) { countsInCol[colsPtr[i]]++; }

  labelVector countsInColOffsets(nCells + 1, other_comm);

  forAll(i, nCells)
  {
    countsInColOffsets[i + 1] = countsInColOffsets[i] + countsInCol[i];
    countsInCol[i] = 0;
  }

  labelVector rowTemp(size, other_comm);
  labelVector colTemp(size, other_comm);

  forAll(i, size)
  {
    label colLocal = colsPtr[i];
    label colStartPos = countsInColOffsets[colLocal];
    label pos = colStartPos + countsInCol[colLocal];

    newOrder[i] = pos;
    rowTemp[pos] = rowsPtr[i];
    colTemp[pos] = colsPtr[i];
    countsInCol[colLocal]++;
  }

  forAll(i, size)
  {
    rowsPtr[i] = rowTemp[i];
    colsPtr[i] = colTemp[i];
  }
}
