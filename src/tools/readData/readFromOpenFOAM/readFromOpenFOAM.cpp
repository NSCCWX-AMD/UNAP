#include "readFromOpenFOAM.hpp"
#include <fstream>
#include <iomanip>
#include <vector>
#include <sstream>

using namespace std;

void UNAP::constructLDUMatrixFromOpenFOAM(lduMatrix& lduA, const char* fileName)
{
	ifstream dataFile;
	dataFile.open(fileName);
    vector<string> vec;
    string temp;

    label nCells = -1;
    label nFaces = -1;
    label nnz = -1;
    int   symm = 0;

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
	        if(pam == 1)
	        {
	            nCells = atoi(s.c_str());
	        }

	        if(pam == 3)
	        {
	            nnz = atoi(s.c_str());
	        }

	        if(pam == 5)
	        {
	            symm = atoi(s.c_str());
	        }
	        pam++;
	    }
    }

    if(symm)
    {
    	nFaces = nnz - nCells;
    }
    else
    {
    	nFaces = (nnz - nCells) / 2;
    }

    if(vec.size() != (nCells + nFaces + 1))
    {
    	UNAPCOUT << "Error in reading " << fileName <<": reading " << vec.size() - 1
    		 << " lines, while nCells = " << nCells
    		 << ", nFaces = " << nFaces << ENDL;

    	ERROR_EXIT;
    }


    labelVector upperAddr(nFaces);
    labelVector lowerAddr(nFaces);
    scalarVector upper(nFaces);
    scalarVector lower(nFaces);
    scalarVector diag(nCells);

    std::vector<std::string>::iterator it;

    int i;
    for(i=0, it=vec.begin(), it++; it!=vec.end(); it++, i++)
    {
    	istringstream is(*it);
	    string s;
	    int pam = 0;
	    if(i < nCells)
	    {
	    	while(is >> s)
		    {
		    	if(pam == 2)
		    	{
		    		diag[i] = atof(s.c_str());
		    	}
		    	pam++;
		    }
	    }
	    else
	    {
	    	label j = i - nCells;
	    	while(is >> s)
		    {
		    	if(pam == 0)
		    	{
		    		lowerAddr[j] = atoi(s.c_str());
		    	}
		    	if(pam == 1)
		    	{
		    		upperAddr[j] = atoi(s.c_str());
		    	}
		    	if(pam == 2)
		    	{
		    		upper[j] = atof(s.c_str());
		    	}
		    	if(!symm)
		    	{
		    		if(pam == 3)
			    	{
			    		lower[j] = atof(s.c_str());
			    	}
		    	}
		    	pam++;
		    }
	    }
    }

    dataFile.close();

    bool reUse = false;

    lduA.setMatrixTopology(upperAddr, lowerAddr, reUse);

    if(symm)
    {
    	lduA.setMatrixCoeffients(diag, upper, upper, reUse);
    }
    else
    {
    	lduA.setMatrixCoeffients(diag, upper, lower, reUse);
    }
}


void UNAP::constructVectorFromOpenFOAM(scalarVector& b, const char* fileName)
{
	ifstream dataFile;
    dataFile.open(fileName);
    vector<string> vec;
    string temp;

    label nCells = -1;

    while(getline(dataFile, temp))
    {
        vec.push_back(temp);
    }

    std::vector<std::string>::iterator it;

    //- read first line
    {
    	it = vec.begin();
	    istringstream is(*it);
	    string s;
	    int pam = 0;
	    while(is >> s)
	    {
	        if(pam == 1)
	        {
	            nCells = atoi(s.c_str());
	        }
	        pam++;
	    }
    }

    if((nCells + 1) != vec.size())
    {
    	UNAPCOUT << "Error in reading " << fileName <<": reading " << vec.size() - 1
    		 << " lines, while nCells = " << nCells << ENDL;
    	ERROR_EXIT;
    }

    if(nCells != b.size())
    {
    	UNAPCOUT << "Error in " << fileName << ": fill size = " << nCells
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
	        if(pam == 0)
	        {
	            b[i] = atof(s.c_str());
	        }
	        pam++;
	    }
    }
    dataFile.close();
}


void UNAP::constructLDUInterfacesFromOpenFOAM(lduMatrix& lduA, const char* fileName)
{
	ifstream dataFile;
    dataFile.open(fileName);
    vector<string> vec;
    string temp;

    while(getline(dataFile, temp))
    {
        vec.push_back(temp);
    }

    label nNeiProcs = -1;
    label nFaces = -1;
    std::vector<std::string>::iterator it;
    //- read first line
    {
    	it = vec.begin();
	    istringstream is(*it);
	    string s;
	    int pam = 0;
	    while(is >> s)
	    {
	        if(pam == 2)
	        {
	            nNeiProcs = atoi(s.c_str());
	        }

	        if(pam == 5)
	        {
	            nFaces = atoi(s.c_str());
	        }
	        pam++;
	    }
    }


    if(vec.size() != (nFaces + nNeiProcs + 1))
    {
    	UNAPCOUT << "Error in reading " << fileName
    		 << ":, reading " << vec.size() - 1
    		 << " lines, while nNeiProcs = " << nNeiProcs
    		 << ", nFaces = " << nFaces << ENDL;

    	ERROR_EXIT;
    }

    label localSize = 0;
    label neighborID = 0;
    PtrList<patch>* patchesPtr = new PtrList<patch>(nNeiProcs);

    forAll(i, nNeiProcs)
    {
    	it++;
    	{
	    	istringstream is(*it);
		    string s;
		    int pam = 0;
		    while(is >> s)
		    {
		        if(pam == 1)
		        {
		            neighborID = atoi(s.c_str());
		        }

		        if(pam == 3)
		        {
		            localSize = atoi(s.c_str());
		        }
		        pam++;
		    }
		}

		patch* patchIPtr = new patch(localSize, MYID, neighborID);
		scalar* localData = new scalar[localSize];
		label*  localFaceCells = new label[localSize];

	    forAll(j, localSize)
	    {
	    	it++;
	    	istringstream is(*it);
		    string s;
		    int pam = 0;

		    while(is >> s)
		    {
		        if(pam == 0)
		        {
		            localFaceCells[j] = atoi(s.c_str());
		        }

		        if(pam == 1)
		        {
		            localData[j] = -atof(s.c_str());
		        }
		        pam++;
		    }
	    }

	    scalarVector* patchCoeffsPtr = new scalarVector(localData, localSize);
   		labelVector* locFaceCellsPtr = new labelVector(localFaceCells, localSize);

   		delete [] localData;
   		delete [] localFaceCells;

   		patchIPtr->patchCoeffs(*patchCoeffsPtr);
   		patchIPtr->faceCells(*locFaceCellsPtr);

   		patchesPtr->setLevel(i, *patchIPtr);
    }

    interfaces* interfacesLocalPtr = new interfaces(*patchesPtr);
	lduA.matrixInterfaces(*interfacesLocalPtr);

    dataFile.close();
}
