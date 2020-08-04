#include "readFromHypre.hpp"

#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

#include "matrixConversion.hpp"
#include "printUNAP.hpp"

using namespace std;

void UNAP::constructLDUMatrixFromHypre(lduMatrix &lduA, const char *fileName)
{
  Communicator *other_comm = lduA.getCommunicator();
  ifstream dataFile;
  dataFile.open(fileName);
  vector<string> vec;
  string temp;

  label rStart = -1, rEnd = -1;
  label cStart = -1, cEnd = -1;

  while (getline(dataFile, temp))
  {
    vec.push_back(temp);
  }

  //- read first line
  {
    vector<string>::iterator it = vec.begin();
    istringstream is(*it);
    string s;
    int pam = 0;
    while (is >> s)
    {
      if (pam == 0)
      {
        rStart = atoi(s.c_str());
      }

      if (pam == 1)
      {
        rEnd = atoi(s.c_str());
      }

      if (pam == 2)
      {
        cStart = atoi(s.c_str());
      }

      if (pam == 3)
      {
        cEnd = atoi(s.c_str());
      }
      pam++;
    }
  }

  //- in fact, there is no need to know end in every processor, except for the
  // last processor
  label nProcs = other_comm->getMySize();
  labelVector procStart(nProcs, other_comm);
  labelVector procEnd(nProcs, other_comm);

  if (nProcs > 1)
  {
    other_comm->barrier();
    other_comm->allGather("allgather1",
                          &rStart,
                          1 * sizeof(label),
                          procStart.values(),
                          1 * sizeof(label));
    other_comm->allGather("allgather2",
                          &rEnd,
                          1 * sizeof(label),
                          procEnd.values(),
                          1 * sizeof(label));
    other_comm->finishTask("allgather1");
    other_comm->finishTask("allgather2");
  }

  const label nCells = rEnd - rStart + 1;

  const label nnzAll = vec.size() - 1;

  labelVector row(nnzAll, other_comm);
  labelVector col(nnzAll, other_comm);
  scalarVector val(nnzAll, other_comm);

  {
    std::vector<std::string>::iterator it;
    int i;
    for (i = 0, it = vec.begin(), it++; it != vec.end(); it++, i++)
    {
      istringstream is(*it);
      string s;
      int pam = 0;
      while (is >> s)
      {
        if (pam == 0)
        {
          row[i] = atoi(s.c_str());
        }

        if (pam == 1)
        {
          col[i] = atoi(s.c_str());
        }

        if (pam == 2)
        {
          val[i] = atof(s.c_str());
        }
        pam++;
      }
    }
  }

  dataFile.close();

  label nnzDiagPart = 0;
  label nnzOffDiagPart = 0;
  label upperSize = 0;

  labelVector procCounts(nProcs, other_comm);

  vector<label> faceToProcNOVec;

  forAll(i, nnzAll)
  {
    if (col[i] >= cStart && col[i] <= cEnd)
    {
      if (row[i] < col[i])
      {
        upperSize++;
      }
      nnzDiagPart++;
    }
    else
    {
      forAll(j, nProcs)
      {
        if (col[i] >= procStart[j] && col[i] <= procEnd[j])
        {
          procCounts[j]++;
          faceToProcNOVec.push_back(j);
        }
      }
      nnzOffDiagPart++;
    }
  }

  label nNeiProcs = 0;
  forAll(i, nProcs)
  {
    if (procCounts[i] > 0)
    {
      nNeiProcs++;
    }
  }

  labelVector nFacesInProc(nNeiProcs, other_comm);
  labelVector neiProcNo(nNeiProcs, other_comm);
  nNeiProcs = 0;
  map<int, int> mapProcID;
  forAll(i, nProcs)
  {
    if (procCounts[i] > 0)
    {
      nFacesInProc[nNeiProcs] = procCounts[i];
      neiProcNo[nNeiProcs] = i;
      mapProcID[i] = nNeiProcs;
      nNeiProcs++;
    }
  }

  labelVector rowLoc(nnzDiagPart, other_comm);
  labelVector colLoc(nnzDiagPart, other_comm);
  scalarVector valLoc(nnzDiagPart, other_comm);

  labelVector rowOffDiag(nnzOffDiagPart, other_comm);
  labelVector colOffDiag(nnzOffDiagPart, other_comm);
  labelVector faceToProcNO(nnzOffDiagPart, other_comm);
  scalarVector valOffDiag(nnzOffDiagPart, other_comm);

  forAll(i, nnzOffDiagPart) { faceToProcNO[i] = faceToProcNOVec[i]; }

  nnzDiagPart = 0;
  nnzOffDiagPart = 0;
  forAll(i, nnzAll)
  {
    if (col[i] >= cStart && col[i] <= cEnd)
    {
      rowLoc[nnzDiagPart] = row[i] - rStart;
      colLoc[nnzDiagPart] = col[i] - cStart;
      valLoc[nnzDiagPart] = val[i];
      nnzDiagPart++;
    }
    else
    {
      rowOffDiag[nnzOffDiagPart] = row[i];
      colOffDiag[nnzOffDiagPart] = col[i];
      valOffDiag[nnzOffDiagPart] = val[i];
      nnzOffDiagPart++;
    }
  }

  reorderCOO(valLoc.values(),
             rowLoc.values(),
             colLoc.values(),
             nCells,
             nnzDiagPart,
             other_comm);

  lduMatrix &diagA = coo2ldu(valLoc.values(),
                             rowLoc.values(),
                             colLoc.values(),
                             nCells,
                             nnzDiagPart,
                             0,
                             other_comm);

  printLDUMatrix(diagA, "A_diag");

  labelVector faceStart(nNeiProcs + 1, other_comm);
  forAll(i, nNeiProcs) { faceStart[i + 1] = faceStart[i] + nFacesInProc[i]; }

  sortInterFaces(valOffDiag,
                 rowOffDiag,
                 colOffDiag,
                 faceToProcNO,
                 nnzOffDiagPart,
                 faceStart,
                 mapProcID,
                 neiProcNo,
                 nNeiProcs,
                 procStart,
                 procEnd);

  bool reUse = false;

  lduA.setMatrixTopology(diagA.upperAddr(), diagA.lowerAddr(), reUse);
  lduA.setMatrixCoeffients(diagA.diag(), diagA.upper(), diagA.lower(), reUse);

  labelVector faceCells(nnzOffDiagPart, other_comm);
  forAll(i, nnzOffDiagPart) { faceCells[i] = rowOffDiag[i] - rStart; }

  if (nProcs > 1)
  {
    constructLDUInterfacesFromHypre(
        lduA, nNeiProcs, neiProcNo, faceStart, faceCells, valOffDiag);
    other_comm->barrier();

    std::cout << "start print interfaces" << ENDL;

    printInterfaces(lduA, "interfaces");
  }
}

void UNAP::sortInterFaces(scalarVector &val,
                          labelVector &row,
                          labelVector &col,
                          const labelVector &faceToProcNO,
                          const label faceSize,
                          const labelVector &faceStart,
                          map<int, int> &mapProcNO,
                          const labelVector &neiProcNo,
                          const label procSize,
                          const labelVector &globalRowStart,
                          const labelVector &globalRowEnd)
{
  Communicator *other_comm = globalRowEnd.getCommunicator();

  label myId = other_comm->getMyId();
  labelVector procCounts(procSize, other_comm);

  scalarVector valTemp(faceSize, other_comm);
  labelVector rowTemp(faceSize, other_comm);
  labelVector colTemp(faceSize, other_comm);

  //- sort as processor
  forAll(i, faceSize)
  {
    label procNO = faceToProcNO[i];
    label procNO_local = mapProcNO[procNO];
    label faceLocation = faceStart[procNO_local] + procCounts[procNO_local];
    valTemp[faceLocation] = val[i];
    rowTemp[faceLocation] = row[i];
    colTemp[faceLocation] = col[i];
    procCounts[procNO_local]++;
  }

  forAll(i, procSize)
  {
    label procNO = neiProcNo[i];
    label localSize = faceStart[i + 1] - faceStart[i];

    scalarVector valTemp2(localSize, other_comm);
    labelVector rowTemp2(localSize, other_comm);
    labelVector colTemp2(localSize, other_comm);

    //- sort lower part
    if (procNO < myId)
    {
      label nCols = globalRowEnd[procNO] - globalRowStart[procNO] + 1;

      labelVector countsInCol(nCols, other_comm);
      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        label colLocal = colTemp[globalPosition] - globalRowStart[procNO];
        countsInCol[colLocal]++;
      }

      labelVector countsInColOffsets(nCols + 1, other_comm);

      forAll(j, nCols)
      {
        countsInColOffsets[j + 1] = countsInColOffsets[j] + countsInCol[j];
        countsInCol[j] = 0;
      }

      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        label colLocal = colTemp[globalPosition] - globalRowStart[procNO];
        label localPosition =
            countsInColOffsets[colLocal] + countsInCol[colLocal];

        valTemp2[localPosition] = valTemp[globalPosition];
        rowTemp2[localPosition] = rowTemp[globalPosition];
        colTemp2[localPosition] = colTemp[globalPosition];
        countsInCol[colLocal]++;
      }

      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        valTemp[globalPosition] = valTemp2[j];
        rowTemp[globalPosition] = rowTemp2[j];
        colTemp[globalPosition] = colTemp2[j];
      }
    }
    //- sort upper part
    else
    {
      label nRows = globalRowEnd[myId] - globalRowStart[myId] + 1;

      labelVector countsInRow(nRows, other_comm);

      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        label rowLocal = rowTemp[globalPosition] - globalRowStart[myId];
        countsInRow[rowLocal]++;
      }

      labelVector countsInRowOffsets(nRows + 1, other_comm);

      forAll(j, nRows)
      {
        countsInRowOffsets[j + 1] = countsInRowOffsets[j] + countsInRow[j];
        countsInRow[j] = 0;
      }

      labelVector posInRow(localSize, other_comm);
      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        label rowLocal = rowTemp[globalPosition] - globalRowStart[myId];

        label localRowSize =
            countsInRowOffsets[rowLocal + 1] - countsInRowOffsets[rowLocal];

        forAll(k, localRowSize)
        {
          if (colTemp[globalPosition] >
              colTemp[faceStart[i] + countsInRowOffsets[rowLocal] + k])
          {
            posInRow[j]++;
          }
        }
      }

      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        label rowLocal = rowTemp[globalPosition] - globalRowStart[myId];
        label localPosition = countsInRowOffsets[rowLocal] + posInRow[j];
        valTemp2[localPosition] = valTemp[globalPosition];
        rowTemp2[localPosition] = rowTemp[globalPosition];
        colTemp2[localPosition] = colTemp[globalPosition];
      }

      forAll(j, localSize)
      {
        label globalPosition = faceStart[i] + j;
        valTemp[globalPosition] = valTemp2[j];
        rowTemp[globalPosition] = rowTemp2[j];
        colTemp[globalPosition] = colTemp2[j];
      }
    }
  }

  forAll(i, faceSize)
  {
    val[i] = valTemp[i];
    row[i] = rowTemp[i];
    col[i] = colTemp[i];
  }
}

void UNAP::constructLDUInterfacesFromHypre(lduMatrix &lduA,
                                           const label nNeiProcs,
                                           const labelVector &destRank,
                                           const labelVector &locPosition,
                                           const labelVector &faceCells,
                                           const scalarVector &data)
{
  Communicator *other_comm = data.getCommunicator();
  label myId = other_comm->getMyId();
  PtrList<patch> *patchesPtr = new PtrList<patch>(nNeiProcs);

  forAll(intI, nNeiProcs)
  {
    const label neighborID = destRank[intI];
    const label localSize = locPosition[intI + 1] - locPosition[intI];
    patch *patchIPtr = new patch(localSize, myId, neighborID);

    scalar *localData = new scalar[localSize];
    label *localFaceCells = new label[localSize];

    forAll(faceI, localSize)
    {
      label start = locPosition[intI] + faceI;
      localData[faceI] = data[start];
      localFaceCells[faceI] = faceCells[start];
    }

    scalarVector *patchCoeffsPtr =
        new scalarVector(localData, localSize, other_comm);
    labelVector *locFaceCellsPtr =
        new labelVector(localFaceCells, localSize, other_comm);

    delete[] localData;
    delete[] localFaceCells;

    patchIPtr->patchCoeffs(*patchCoeffsPtr);
    patchIPtr->faceCells(*locFaceCellsPtr);

    patchesPtr->setLevel(intI, *patchIPtr);
  }

  interfaces *interfacesLocalPtr = new interfaces(*patchesPtr, other_comm);
  lduA.matrixInterfaces(*interfacesLocalPtr);
}

void UNAP::constructVectorFromHypre(scalarVector &b, const char *fileName)
{
  ifstream dataFile;
  dataFile.open(fileName);
  vector<string> vec;
  string temp;

  while (getline(dataFile, temp))
  {
    vec.push_back(temp);
  }

  std::vector<std::string>::iterator it;

  label rStart = -1, rEnd = -1;

  //- read first line
  {
    it = vec.begin();
    istringstream is(*it);
    string s;
    int pam = 0;
    while (is >> s)
    {
      if (pam == 0)
      {
        rStart = atoi(s.c_str());
      }

      if (pam == 1)
      {
        rEnd = atoi(s.c_str());
      }
      pam++;
    }
  }

  label nCells = rEnd - rStart + 1;

  if ((nCells + 1) != vec.size())
  {
    std::cout << "Error in reading " << fileName << ": reading "
              << vec.size() - 1 << " lines, while nCells = " << nCells << ENDL;
    ERROR_EXIT;
  }

  if (nCells != b.size())
  {
    std::cout << "Error in " << fileName << ": fill size = " << nCells
              << ", while allocated size = " << b.size() << ENDL;
    ERROR_EXIT;
  }

  int i;
  for (i = 0, it = vec.begin(), it++; it != vec.end(); it++, i++)
  {
    istringstream is(*it);
    string s;
    int pam = 0;
    while (is >> s)
    {
      if (pam == 1)
      {
        b[i] = atof(s.c_str());
      }
      pam++;
    }
  }
  dataFile.close();
}
