/*
	author: Jin yujun
	datatime: 2019/8/19
	email: miluan@bupt.cdu.cn
*/
#include <bits/stdc++.h>
#include <iostream>
#include "cooRCM.hpp"
#include "rcm.h"

using namespace COORCM;

#define DELETE_POINTER(ptr) if(ptr) \
{ \
    delete [] ptr; \
    ptr = NULL; \
}

#define DELETE_OBJECT_POINTER(ptr) if(ptr) \
{ \
    delete ptr; \
    ptr = NULL; \
}

int findIndex(std::vector<std::pair<int, int> > a, int x)
{
    int nRet = -1;
	for (int i = 0; i < a.size(); i++)
	{
		if (a[i].first == x)
		{
            nRet = i;
            break;
        }
    }
	return nRet;
}


ReorderingSSM::ReorderingSSM(const int nnz, const int nRows, const int* row, const int* col)
:
	_nnz(nnz),
	_nRows(nRows),
	_row(row),
	_col(col)
{
	_vetexOrderPtr = new std::vector<int>();

	_edgeOrderPtr = new std::vector<int>();

	int Max = 0;
	for(int i=0; i<nnz ; i++)
	{
		int temp = abs(row[i] - col[i]);
		if(temp > Max)
		{
			Max = temp;
		}
	}
	std::cout << "Bandwidth before reordering: " << Max << std::endl;

	ReverseCuthillMckee();
}

ReorderingSSM::~ReorderingSSM()
{
	DELETE_OBJECT_POINTER(_vetexOrderPtr);
	DELETE_OBJECT_POINTER(_edgeOrderPtr);
}

std::vector<int> ReorderingSSM::_globalDegree;

bool ReorderingSSM::compareDegree(int i, int j)
{
	return _globalDegree[i] - _globalDegree[j];
}

void ReorderingSSM::degreeGenerator()
{
	std::cout << "in degreeGenerator" << std::endl;
	_globalDegree.clear();
	int* nColsInRow = new int[_nRows];
	for(int i=0; i<_nRows; i++)
	{
		nColsInRow[i] = 0;
	}
	for(int i=0; i<_nnz; i++)
	{
		nColsInRow[_row[i]]++;
	}
	for(int i=0; i<_nRows; i++)
	{
		_globalDegree.push_back(nColsInRow[i]);
	}
	DELETE_POINTER(nColsInRow);

	std::cout << "end degreeGenerator" << std::endl;
}

void ReorderingSSM::CuthillMckee()
{
	std::cout << "in CuthillMckee" << std::endl;
	degreeGenerator();
	std::queue<int>* Qptr = new std::queue<int>();
	std::queue<int>& Q = *Qptr;

	std::vector<std::pair<int, int> >* notVisitedPtr = new std::vector<std::pair<int, int> >();
	std::vector<std::pair<int, int> >& notVisited = *notVisitedPtr;

	std::vector<int>& R = *_vetexOrderPtr;

	for (int i = 0; i < _globalDegree.size(); i++)
	{
		notVisited.push_back(std::make_pair(i, _globalDegree[i]));
	}

	// Vector notVisited helps in running BFS
	// even when there are dijoind graphs
	while (!notVisited.empty())
	{
		int minNodeIndex = 0;

		for (int i = 0; i < notVisited.size(); i++)
		{
			if (notVisited[i].second < notVisited[minNodeIndex].second)
			{
				minNodeIndex = i;
			}
		}

		Q.push(notVisited[minNodeIndex].first);
		notVisited.erase(notVisited.begin() +// minNodeIndex)
										 findIndex(notVisited,
											notVisited[Q.front()].first));
		// Simple BFS
		while (!Q.empty())
		{
			std::vector<int>* toSortPtr = new std::vector<int>();
			std::vector<int>& toSort = *toSortPtr;
			int y = Q.front();

			for(int i=0; i<_nnz; i++)
			{
				if(_row[i]==Q.front() && _row[i]!=_col[i] && findIndex(notVisited, _col[i]) != -1)
				{
					toSort.push_back(_col[i]);
					notVisited.erase(notVisited.begin() + findIndex(notVisited, _col[i]));
				}
			}
			sort(toSort.begin(), toSort.end(), compareDegree);

			for (int i = 0; i < toSort.size(); i++)
			{
				Q.push(toSort[i]);
			}

			R.push_back(Q.front());
			Q.pop();
		}
	}

	DELETE_OBJECT_POINTER(Qptr);
	DELETE_OBJECT_POINTER(notVisitedPtr);
	std::cout << "end CuthillMckee" << std::endl;
}


void ReorderingSSM::ReverseCuthillMckee()
{
	std::cout << "in ReverseCuthillMckee" << std::endl;
	CuthillMckee();
	std::vector<int>& cuthill = *_vetexOrderPtr;

	int n = cuthill.size();

	if (n % 2 == 0)
	{
		n -= 1;
	}

	n = n / 2;

	//- reverse
	for (int i = 0; i <= n; i++)
	{
		int j = cuthill[cuthill.size() - 1 - i];
		cuthill[cuthill.size() - 1 - i] = cuthill[i];
		cuthill[i] = j;
	}

	//-
	int* ord = new int[_nRows];
	for(int i=0; i<_nRows; i++)
	{
		ord[i] = i;
	}
	int otto;
	for (int i = 0; i < _nRows - 1; i++)
	{
        for (int j = _nRows - 1; j > i; j--)
		{
        	if (cuthill[j - 1] > cuthill[j])
        	{
       			otto = cuthill[j - 1];
       			cuthill[j - 1] = cuthill[j];
       			cuthill[j] = otto;

				otto = ord[j-1];
				ord[j-1] = ord[j];
				ord[j] = otto;
            }
        }
    }
	for(int i=0; i<_nRows; i++)
	{
		cuthill[i] = ord[i];
	}

	int Max = 0;
	for(int i=0; i<_nnz; i++)
   	{
   	 	int newrow = cuthill[_row[i]];
   	 	int newcol = cuthill[_col[i]];

   	 	int otto = abs(newrow - newcol);
        if(otto > Max)
        {
            Max = otto;
        }
   	}
   	std::cout << "Bandwidth after  reordering: " << Max << std::endl;
	DELETE_POINTER(ord);
	std::cout << "end ReverseCuthillMckee" << std::endl;
}

std::vector<int>* ReorderingSSM::getVetexOrder()
{
	return _vetexOrderPtr;
}

std::vector<int>* ReorderingSSM::getEdgeOrder()
{
	//need to change
	std::vector<int>& edgeorder = *_edgeOrderPtr;

	int otto;
	int* newrow = new int[_nnz];
	int* newcol = new int[_nnz];
	int* order = new int[_nnz];

	for(int i=0; i<_nnz; i++)
	{
		edgeorder.push_back(i);
		order[i] = i;
		newrow[i] = (*_vetexOrderPtr)[_row[i]];
   	 	newcol[i] = (*_vetexOrderPtr)[_col[i]];
	}




	for (int i = 0; i < _nnz - 1; i++)
	{
        for (int j = _nnz - 1; j > i; j--)
		{
        	if (newcol[j - 1] > newcol[j])
        	{
       			otto = newcol[j - 1];
       			newcol[j - 1] = newcol[j];
       			newcol[j] = otto;

				otto = newrow[j-1];
				newrow[j-1] = newrow[j];
				newrow[j] = otto;

				otto = order[j-1];
				order[j-1] = order[j];
				order[j] = otto;
            }
        }
    }
	for (int i = 0; i < _nnz - 1; i++)
    {
        for (int j = _nnz - 1; j > i; j--)
        {
            if (newrow[j - 1] > newrow[j])
            {
	            otto = newcol[j - 1];
	            newcol[j - 1] = newcol[j];
	            newcol[j] = otto;
	            otto = newrow[j-1];
	            newrow[j-1] = newrow[j];
	            newrow[j] = otto;

	            otto = order[j-1];
				order[j-1] = order[j];
				order[j] = otto;
            }
        }
    }
    for (int i = 0; i < _nnz - 1; i++)
    {
        for (int j = _nnz - 1; j > i; j--)
        {
            if (order[j - 1] > order[j])
            {
               	otto = order[j-1];
               	order[j-1] = order[j];
               	order[j] = otto;

            	otto = edgeorder[j-1];
				edgeorder[j-1] = edgeorder[j];
				edgeorder[j] = otto;
            }
        }
    }
    DELETE_POINTER(newrow);
	DELETE_POINTER(newcol);
	DELETE_POINTER(order);
	return _edgeOrderPtr;
}

void rcmCOO_nowrite(const int nnz,
					const int nRows,
					const int *row,
					const int *col,
					int *postVetexOrder,
					int *postValOrder)
{
	ReorderingSSM m(nnz, nRows, row, col);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);
	std::vector<int>* EdgeorderPtr = m.getEdgeOrder();
	std::vector<int>& Edgeorder = *(EdgeorderPtr);

	for(int ij=0; ij < Cellorder.size(); ij++)
	{
		postVetexOrder[ij] = Cellorder[ij];
	}
	for(int ik=0; ik < Edgeorder.size(); ik++)
	{
		postValOrder[ik] = Edgeorder[ik];
	}
}

void rcmCOO_rewrite(const int nnz, const int nRows, int *row, int *col, double *val)
{
	ReorderingSSM m(nnz, nRows, row, col);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);
	std::vector<int>* EdgeorderPtr = m.getEdgeOrder();
	std::vector<int>& Edgeorder = *(EdgeorderPtr);

	int* ncol = new int[nnz];
	int* nrow = new int[nnz];
	double* newval = new double[nnz];
    for(int j=0; j<nnz; j++)
    {
    	nrow[Edgeorder[j]] = Cellorder[row[j]];
    	ncol[Edgeorder[j]] = Cellorder[col[j]];
        newval[Edgeorder[j]] = val[j];
    }
    for(int j=0; j<nnz; j++)
    {
        col[j] = ncol[j];
        row[j] = nrow[j];
		val[j] = newval[j];
    }

	DELETE_POINTER(ncol);
	DELETE_POINTER(nrow);
    DELETE_POINTER(newval);
}

void rcmCSR_nowrite(const int nnz,
			 		const int nRows,
			 		const int* rowsOffset,
			 		const int* col,
			 		int* postVetexOrder,
			 		int* postValOrder,
			 		int* newRowsOffset)
{
	int* row = new int[nnz];
	int j = 0;
	for(int i=0; i<nnz; i++)
	{
		if(i >= rowsOffset[j+1])
		{
			j++;
		}
		row[i] = j;
	}

	ReorderingSSM m(nnz, nRows, row, col);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);
	std::vector<int>* EdgeorderPtr = m.getEdgeOrder();
	std::vector<int>& Edgeorder = *(EdgeorderPtr);

	for(int ij=0; ij < Cellorder.size(); ij++)
	{
		postVetexOrder[ij] = Cellorder[ij];
	}
	for(int ik=0; ik < Edgeorder.size(); ik++)
	{
		postValOrder[ik] = Edgeorder[ik];
	}
	newRowsOffset[0] = 0;

	for(int j=1; j<nRows+1; j++)
	{
		newRowsOffset[Cellorder[j-1]+1] = rowsOffset[j] - rowsOffset[j-1];
	}
	for(int j=1; j<nRows+1; j++)
	{
		newRowsOffset[j] += newRowsOffset[j-1];
	}

	DELETE_POINTER(row);
}

void rcmCSR_rewrite(const int nnz, const int nRows, int* rowsOffset, int* col, double* value)
{
	int* row = new int[nnz];
	int k = 0;
	for(int i=0; i<nnz; i++)
	{
		if(i >= rowsOffset[k+1])
		{
			k++;
		}
		row[i] = k;
	}
	ReorderingSSM m(nnz, nRows, row, col);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);
	std::vector<int>* EdgeorderPtr = m.getEdgeOrder();
	std::vector<int>& Edgeorder = *(EdgeorderPtr);

	int* newcol = new int[nnz];
	double* newval = new double[nnz];

	for(int j=0; j<nnz; j++)
	{
		newcol[Edgeorder[j]] = Cellorder[col[j]];
		newval[Edgeorder[j]] = value[j];
	}

	for(int j=0; j<nnz; j++)
	{
		col[j] = newcol[j];
		value[j] = newval[j];
	}

	int* newRowsOffset = new int[nRows+1];
	newRowsOffset[0] = 0;
	for(int j=1; j<nRows+1; j++)
	{
		newRowsOffset[Cellorder[j-1]+1] = rowsOffset[j] - rowsOffset[j-1];
	}
	for(int j=1; j<nRows+1; j++)
	{
		rowsOffset[j] = newRowsOffset[j] + newRowsOffset[j-1];
	}

	DELETE_POINTER(row);
	DELETE_POINTER(newRowsOffset);
    DELETE_POINTER(newcol);
    DELETE_POINTER(newval);
}

void rcmLDU_nowrite(const int nnz,
					const int nRows,
					const int* row,
					const int* col,
					int* postVetexOrder,
					int* postEdgeOrder)
{
	int len = 2*nnz;
	int* col2 = new int[len];
	int* row2 = new int[len];
	int* nColsInRow = new int[nRows];
	int* colsOffset = new int[nRows+1];
	colsOffset[0] = 0;
	for(int i=0; i<nRows; i++)
	{
		nColsInRow[i] = 0;
		colsOffset[i+1] = 0;
	}
	for(int i=0; i<nnz; i++)
	{
		int r=row[i];
		int c=col[i];

		nColsInRow[r]++;
		nColsInRow[c]++;
	}

	for(int i=0; i<nRows; i++)
	{
		colsOffset[i+1] = colsOffset[i] + nColsInRow[i];
		nColsInRow[i] = 0;
	}

	for(int i=0; i<nnz; i++)
	{
		int r=row[i];
		int c=col[i];

		row2[colsOffset[r]+nColsInRow[r]] = r;
		col2[colsOffset[r]+nColsInRow[r]] = c;
		nColsInRow[r]++;

		row2[colsOffset[c]+nColsInRow[c]] = c;
		col2[colsOffset[c]+nColsInRow[c]] = r;
		nColsInRow[c]++;
	}
	DELETE_POINTER(nColsInRow);
	DELETE_POINTER(colsOffset);

	ReorderingSSM m(len, nRows, row2, col2);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);
	int* order = new int[nnz];

	for(int it=0; it<nRows; it++)
	{
		postVetexOrder[it] = Cellorder[it];
	}

	for(int it=0; it<nnz; it++)
	{
		postEdgeOrder[it] = it;
		order[it] = it;
	}

	//need to change
	int* nrow = new int[nnz];
	int* ncol = new int[nnz];

	int otto;
	for(int i=0; i<nnz; i++)
    {
    	nrow[i] = Cellorder[row[i]];
    	ncol[i] = Cellorder[col[i]];
    	if(nrow[i] > ncol[i])
    	{
    		order[i] = -order[i]-1;
    		otto = nrow[i];
			nrow[i] = ncol[i];
			ncol[i] = otto;
		}
    }

    for(int i=0; i<nnz; i++)
    {

    }


	for (int i = 0; i < nnz - 1; i++)
	{
        for (int j = nnz - 1; j > i; j--)
		{
        	if (ncol[j - 1] > ncol[j])
        	{
       			otto = ncol[j - 1];
       			ncol[j-1] = ncol[j];
       			ncol[j] = otto;

				otto = nrow[j-1];
				nrow[j-1] = nrow[j];
				nrow[j] = otto;

				otto = order[j-1];
				order[j-1] = order[j];
				order[j] = otto;
            }
        }
    }
	for (int i = 0; i < nnz - 1; i++)
    {
        for (int j = nnz - 1; j > i; j--)
        {
            if (nrow[j - 1] > nrow[j])
            {
               otto = ncol[j - 1];
               ncol[j - 1] = ncol[j];
               ncol[j] = otto;

            	otto = nrow[j-1] ;
                nrow[j-1] = nrow[j];
                nrow[j] = otto;

               	otto = order[j-1];
				order[j-1] = order[j];
				order[j] = otto;
            }
         }
    }

    for (int i = 0; i < nnz ; i++)
    {
    	if(order[i] < 0)
    	{
    		order[i] = -order[i] - 1;
    		postEdgeOrder[i] = -postEdgeOrder[i] -1;
		}
	}
    for (int i = 0; i < nnz - 1; i++)
    {
        for (int j = nnz - 1; j > i; j--)
        {
            if (order[j - 1] > order[j])
            {
               otto = order[j - 1];
               order[j - 1] = order[j];
               order[j] = otto;

               	otto = postEdgeOrder[j-1];
				postEdgeOrder[j-1] = postEdgeOrder[j];
				postEdgeOrder[j] = otto;
            }
         }
    }

	DELETE_POINTER(ncol);
	DELETE_POINTER(nrow);
	DELETE_POINTER(col2);
	DELETE_POINTER(row2);
}

void rcmLDU_rewrite(const int nnz,
					const int nRows,
					int* row,
					int* col,
					double* diagVal,
					double* upperVal,
					double* lowerVal)
{
	int len = 2*nnz;
	int* col2 = new int[len];
	int* row2 = new int[len];
	int* nColsInRow = new int[nRows];
	int* colsOffset = new int[nRows+1];
	colsOffset[0] = 0;
	for(int i=0; i<nRows; i++)
	{
		nColsInRow[i] = 0;
		colsOffset[i+1] = 0;
	}
	for(int i=0; i<nnz; i++)
	{
		int r=row[i];
		int c=col[i];

		nColsInRow[r]++;
		nColsInRow[c]++;
	}

	for(int i=0; i<nRows; i++)
	{
		colsOffset[i+1] = colsOffset[i] + nColsInRow[i];
		nColsInRow[i] = 0;
	}

	for(int i=0; i<nnz; i++)
	{
		int r=row[i];
		int c=col[i];

		row2[colsOffset[r]+nColsInRow[r]] = r;
		col2[colsOffset[r]+nColsInRow[r]] = c;
		nColsInRow[r]++;

		row2[colsOffset[c]+nColsInRow[c]] = c;
		col2[colsOffset[c]+nColsInRow[c]] = r;
		nColsInRow[c]++;
	}
	DELETE_POINTER(nColsInRow);
	DELETE_POINTER(colsOffset);

	ReorderingSSM m(len, nRows, row2, col2);

	std::vector<int>* CellorderPtr = m.getVetexOrder();
	std::vector<int>& Cellorder = *(CellorderPtr);

	int temp;
	double tem;
	int* newcol = new int[nnz];
	int* newrow = new int[nnz];
	double* newdiag = new double[nRows];

	for(int i=0; i<nRows; i++)
	{
		newdiag[Cellorder[i]] = diagVal[i];
	}
	for(int i=0;i<nRows;i++)
	{
		diagVal[i] = newdiag[i];
	}
	for(int i=0; i<nnz; i++)
    {
    	newrow[i] = Cellorder[row[i]];
    	newcol[i] = Cellorder[col[i]];
    }

	for(int j=0; j<nnz; j++)
	{
		if(newrow[j] > newcol[j])
		{
			tem = upperVal[j];
			upperVal[j] = lowerVal[j];
			lowerVal[j] = tem;
			temp = newrow[j];
			newrow[j] = newcol[j];
			newcol[j] = temp;
		}
	}
	for(int j=0; j<nnz; j++)
    {
 		col[j] = newcol[j];
 		row[j] = newrow[j];
 	}

	int otto;
	double ott;
	for (int i = 0; i < nnz - 1; i++)
	{
        for (int j = nnz - 1; j > i; j--)
		{
        	if (col[j - 1] > col[j])
        	{
       			otto = col[j - 1];
       			col[j - 1] = col[j];
       			col[j] = otto;

				otto = row[j-1];
				row[j-1] = row[j];
				row[j] = otto;

				ott = upperVal[j-1];
				upperVal[j-1] = upperVal[j];
				upperVal[j] = ott;
				ott = lowerVal[j-1];
				lowerVal[j-1] = lowerVal[j];
				lowerVal[j] = ott;
            }
        }
    }
	for (int i = 0; i < nnz - 1; i++)
    {
        for (int j = nnz - 1; j > i; j--)
        {
            if (row[j - 1] > row[j])
            {
               otto = col[j - 1];
               col[j - 1] = col[j];
               col[j] = otto;
               otto = row[j-1];
               row[j-1] = row[j];
               row[j] = otto;
               ott = upperVal[j-1];
               upperVal[j-1] = upperVal[j];
               upperVal[j] = ott;
               ott = lowerVal[j-1];
               lowerVal[j-1] = lowerVal[j];
               lowerVal[j] = ott;
            }
         }
    }

	DELETE_POINTER(newdiag);
    DELETE_POINTER(col2);
	DELETE_POINTER(row2);
	DELETE_POINTER(newcol);
	DELETE_POINTER(newrow);
}

