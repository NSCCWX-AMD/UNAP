#include "readFromHypre.hpp"
#include "matrixConversion.hpp"
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>
#include "printUNAP.hpp"

using namespace std;

void UNAP::constructLDUMatrixFromHypre(lduMatrix& lduA, const char* fileName)
{
	ifstream dataFile;
	dataFile.open(fileName);
    vector<string> vec;
    string temp;

    label rStart = -1, rEnd = -1;
    label cStart = -1, cEnd = -1;

    while(getline(dataFile, temp))
    {
        vec.push_back(temp);
    }

    //- read first line
    {
    	vector<string>::iterator it = vec.begin();
	    istringstream is(*it);
	    string s;
	    int pam = 0;
	    while(is >> s)
	    {
	        if(pam == 0)
	        {
	            rStart = atoi(s.c_str());
	        }

	        if(pam == 1)
	        {
	            rEnd = atoi(s.c_str());
	        }

	        if(pam == 2)
	        {
	            cStart = atoi(s.c_str());
	        }

	        if(pam == 3)
	        {
	            cEnd = atoi(s.c_str());
	        }
	        pam++;
	    }
    }

    //- in fact, there is no need to know end in every processor, except for the last processor
    labelField procStart(NPROCS);
    labelField procEnd(NPROCS);

    if(PARRUN)
    {
    	MPI_Barrier(MPI_COMM_WORLD);
    	MPI_Allgather(&rStart, 1, MPI_LABEL, procStart.values(), 1, MPI_LABEL, MPI_COMM_WORLD);
    	MPI_Allgather(&rEnd, 1, MPI_LABEL, procEnd.values(), 1, MPI_LABEL, MPI_COMM_WORLD);
    }

    const label nCells = rEnd - rStart + 1;

    const label nnzAll = vec.size() - 1;

    labelField  row(nnzAll);
    labelField  col(nnzAll);
    scalarField val(nnzAll);

    {
    	std::vector<std::string>::iterator it;
	    int i;
	    for(i=0, it=vec.begin(), it++; it!=vec.end(); it++, i++)
	    {
	    	istringstream is(*it);
		    string s;
		    int pam = 0;
	    	while(is >> s)
		    {
		    	if(pam == 0)
		    	{
		    		row[i] = atoi(s.c_str());
		    	}

		    	if(pam == 1)
		    	{
		    		col[i] = atoi(s.c_str());
		    	}

		    	if(pam == 2)
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

    labelField procCounts(NPROCS);

    vector<label> faceToProcNOVec;

    forAll(i, nnzAll)
    {
    	if(col[i]>=cStart && col[i]<=cEnd)
    	{
    		if(row[i] < col[i])
    		{
    			upperSize++;
    		}
    		nnzDiagPart++;
    	}
    	else
    	{
    		forAll(j, NPROCS)
    		{
    			if(col[i]>=procStart[j] && col[i]<=procEnd[j])
    			{
    				procCounts[j]++;
    				faceToProcNOVec.push_back(j);
    			}
    		}
    		nnzOffDiagPart++;
    	}
    }

    label nNeiProcs = 0;
    forAll(i, NPROCS)
    {
    	if(procCounts[i] > 0)
    	{
    		nNeiProcs++;
    	}
    }

    labelField nFacesInProc(nNeiProcs);
    labelField neiProcNo(nNeiProcs);
    nNeiProcs = 0;
    map<int, int> mapProcID;
    forAll(i, NPROCS)
    {
    	if(procCounts[i] > 0)
    	{
    		nFacesInProc[nNeiProcs] = procCounts[i];
    		neiProcNo[nNeiProcs] = i;
    		mapProcID[i] = nNeiProcs;
    		nNeiProcs++;
    	}
    }


    labelField  rowLoc(nnzDiagPart);
    labelField  colLoc(nnzDiagPart);
    scalarField valLoc(nnzDiagPart);

    labelField  rowOffDiag(nnzOffDiagPart);
    labelField  colOffDiag (nnzOffDiagPart);
   	labelField  faceToProcNO(nnzOffDiagPart);
   	scalarField valOffDiag(nnzOffDiagPart);

   	forAll(i, nnzOffDiagPart)
    {
    	faceToProcNO[i] = faceToProcNOVec[i];
    }

    nnzDiagPart = 0;
    nnzOffDiagPart = 0;
    forAll(i, nnzAll)
    {
    	if(col[i]>=cStart && col[i]<=cEnd)
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

    reorderCOO(valLoc.values(), rowLoc.values(), colLoc.values(), nCells, nnzDiagPart);

    lduMatrix& diagA = coo2ldu(
    							valLoc.values(),
    							rowLoc.values(),
    							colLoc.values(),
    							nCells,
    							nnzDiagPart,
    							0
    						   );

    printLDUMatrix(diagA, "A_diag");

    labelField faceStart(nNeiProcs+1);
	forAll(i, nNeiProcs)
	{
		faceStart[i+1] = faceStart[i] + nFacesInProc[i];
	}

    sortInterFaces
    (
    	valOffDiag,
    	rowOffDiag,
    	colOffDiag,
    	faceToProcNO,
    	nnzOffDiagPart,
    	faceStart,
    	mapProcID,
    	neiProcNo,
    	nNeiProcs,
    	procStart,
    	procEnd
    );

    bool reUse = false;

    lduA.setMatrixTopology(diagA.upperAddr(), diagA.lowerAddr(), reUse);
    lduA.setMatrixCoeffients(diagA.diag(), diagA.upper(), diagA.lower(), reUse);

    labelField faceCells(nnzOffDiagPart);
    forAll(i, nnzOffDiagPart)
    {
    	faceCells[i] = rowOffDiag[i] - rStart;
    }

    if(PARRUN)
    {
    	constructLDUInterfacesFromHypre(lduA, nNeiProcs, neiProcNo, faceStart, faceCells, valOffDiag);
    	MPI_Barrier(MPI_COMM_WORLD);
    	if(!MYID)
    	{
    		COUT << "start print interfaces" << ENDL;
    	}
    	printInterfaces(lduA, "interfaces");
    }
}


void UNAP::sortInterFaces
(
	scalarField& val,
	labelField&  row,
	labelField&  col,
	const labelField& faceToProcNO,
	const label  faceSize,
	const labelField& faceStart,
	map<int, int>& mapProcNO,
	const labelField& neiProcNo,
	const label  procSize,
	const labelField& globalRowStart,
	const labelField& globalRowEnd
)
{
	labelField procCounts(procSize);

	scalarField valTemp(faceSize);
	labelField  rowTemp(faceSize);
	labelField  colTemp(faceSize);

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
		label localSize = faceStart[i+1] - faceStart[i];

		scalarField valTemp2(localSize);
		labelField  rowTemp2(localSize);
		labelField  colTemp2(localSize);

		//- sort lower part
		if(procNO < MYID)
		{
			label nCols = globalRowEnd[procNO] - globalRowStart[procNO] + 1;

			labelField  countsInCol(nCols);
			forAll(j, localSize)
			{
				label globalPosition = faceStart[i] + j;
				label colLocal = colTemp[globalPosition] - globalRowStart[procNO];
				countsInCol[colLocal]++;
			}

			labelField countsInColOffsets(nCols+1);

			forAll(j, nCols)
			{
				countsInColOffsets[j+1] = countsInColOffsets[j] + countsInCol[j];
				countsInCol[j] = 0;
			}

			forAll(j, localSize)
			{
				label globalPosition = faceStart[i] + j;
				label colLocal = colTemp[globalPosition] - globalRowStart[procNO];
				label localPosition = countsInColOffsets[colLocal] + countsInCol[colLocal];

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
			label nRows = globalRowEnd[MYID] - globalRowStart[MYID] + 1;

			labelField  countsInRow(nRows);

			forAll(j, localSize)
			{
				label globalPosition = faceStart[i] + j;
				label rowLocal = rowTemp[globalPosition] - globalRowStart[MYID];
				countsInRow[rowLocal]++;
			}


			labelField countsInRowOffsets(nRows+1);

			forAll(j, nRows)
			{
				countsInRowOffsets[j+1] = countsInRowOffsets[j] + countsInRow[j];
				countsInRow[j] = 0;
			}


			labelField posInRow(localSize);
			forAll(j, localSize)
			{
				label globalPosition = faceStart[i] + j;
				label rowLocal = rowTemp[globalPosition] - globalRowStart[MYID];

				label localRowSize = countsInRowOffsets[rowLocal+1] - countsInRowOffsets[rowLocal];

				forAll(k, localRowSize)
				{
					if(colTemp[globalPosition] > colTemp[faceStart[i]+countsInRowOffsets[rowLocal]+k])
					{
						posInRow[j]++;
					}
				}
			}

			forAll(j, localSize)
			{
				label globalPosition = faceStart[i] + j;
				label rowLocal = rowTemp[globalPosition] - globalRowStart[MYID];
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


void UNAP::constructLDUInterfacesFromHypre
(
	lduMatrix&  lduA,
	const label nNeiProcs,
	const labelField&  destRank,
	const labelField&  locPosition,
	const labelField&  faceCells,
	const scalarField& data
)
{
	PtrList<patch>* patchesPtr = new PtrList<patch>(nNeiProcs);

	forAll(intI, nNeiProcs)
	{
		const label neighborID = destRank[intI];
		const label localSize = locPosition[intI+1] - locPosition[intI];
		patch* patchIPtr = new patch(localSize, MYID, neighborID);

		scalar* localData = new scalar[localSize];
		label*  localFaceCells = new label[localSize];

		forAll(faceI, localSize)
		{
			label start = locPosition[intI] + faceI;
			localData[faceI] = data[start];
			localFaceCells[faceI] = faceCells[start];
		}

		scalarField* patchCoeffsPtr = new scalarField(localData, localSize);
   		labelField* locFaceCellsPtr = new labelField(localFaceCells, localSize);

   		delete [] localData;
   		delete [] localFaceCells;

   		patchIPtr->patchCoeffs(*patchCoeffsPtr);
   		patchIPtr->faceCells(*locFaceCellsPtr);

   		patchesPtr->setLevel(intI, *patchIPtr);
	}

	interfaces* interfacesLocalPtr = new interfaces(*patchesPtr);
	lduA.matrixInterfaces(*interfacesLocalPtr);
}


void UNAP::constructVectorFromHypre(scalarField& b, const char* fileName)
{
	ifstream dataFile;
    dataFile.open(fileName);
    vector<string> vec;
    string temp;

    while(getline(dataFile, temp))
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
	    while(is >> s)
	    {
	        if(pam == 0)
	        {
	            rStart = atoi(s.c_str());
	        }

	        if(pam == 1)
	        {
	            rEnd = atoi(s.c_str());
	        }
	        pam++;
	    }
    }

    label nCells = rEnd - rStart + 1;

    if((nCells + 1) != vec.size())
    {
    	COUT << "Error in reading " << fileName <<": reading " << vec.size() - 1
    		 << " lines, while nCells = " << nCells << ENDL;
    	ERROR_EXIT;
    }

    if(nCells != b.size())
    {
    	COUT << "Error in " << fileName << ": fill size = " << nCells
    		 << ", while allocated size = " << b.size() << ENDL;
    	ERROR_EXIT;
    }

    int i;
    for(i=0, it=vec.begin(), it++; it!=vec.end(); it++, i++)
    {
    	istringstream is(*it);
	    string s;
	    int pam = 0;
	    while(is >> s)
	    {
	        if(pam == 1)
	        {
	            b[i] = atof(s.c_str());
	        }
	        pam++;
	    }
    }
    dataFile.close();
}
