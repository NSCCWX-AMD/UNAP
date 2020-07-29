/*
  author: Jin yujun
  datatime: 2019/8/19
  email: miluan@bupt.cdu.cn
*/
#include "cooRCM.hpp"

#include <algorithm>
#include <iostream>
#include <queue>

#include "rcmf.h"

using namespace COORCM;

label findIndex(std::vector<std::pair<label, label> > a, label x)
{
  label nRet = -1;
  for (label i = 0; i < a.size(); i++)
  {
    if (a[i].first == x)
    {
      nRet = i;
      break;
    }
  }
  return nRet;
}

ReorderingSSM::ReorderingSSM(const label nnz,
                             const label nRows,
                             const label *row,
                             const label *col)
    : _nnz(nnz), _nRows(nRows), _row(row), _col(col)
{
  _vetexOrderPtr = new std::vector<label>();

  _edgeOrderPtr = new std::vector<label>();

  label Max = 0;
  for (label i = 0; i < nnz; i++)
  {
    label temp = abs(row[i] - col[i]);
    if (temp > Max)
    {
      Max = temp;
    }
  }
  POUT << "Bandwidth before reordering: " << Max << std::endl;

  ReverseCuthillMckee();
}

ReorderingSSM::~ReorderingSSM()
{
  DELETE_OBJECT_POINTER(_vetexOrderPtr);
  DELETE_OBJECT_POINTER(_edgeOrderPtr);
}

std::vector<label> ReorderingSSM::_globalDegree;

bool ReorderingSSM::compareDegree(label i, label j)
{
  return _globalDegree[i] - _globalDegree[j];
}

void ReorderingSSM::degreeGenerator()
{
  POUT << "in degreeGenerator" << std::endl;
  _globalDegree.clear();
  label *nColsInRow = new label[_nRows];
  for (label i = 0; i < _nRows; i++)
  {
    nColsInRow[i] = 0;
  }
  for (label i = 0; i < _nnz; i++)
  {
    nColsInRow[_row[i]]++;
  }
  for (label i = 0; i < _nRows; i++)
  {
    _globalDegree.push_back(nColsInRow[i]);
  }
  DELETE_POINTER(nColsInRow);

  POUT << "end degreeGenerator" << std::endl;
}

void ReorderingSSM::CuthillMckee()
{
  POUT << "in CuthillMckee" << std::endl;
  degreeGenerator();
  std::queue<label> *Qptr = new std::queue<label>();
  std::queue<label> &Q = *Qptr;

  std::vector<std::pair<label, label> > *notVisitedPtr =
      new std::vector<std::pair<label, label> >();
  std::vector<std::pair<label, label> > &notVisited = *notVisitedPtr;

  std::vector<label> &R = *_vetexOrderPtr;

  for (label i = 0; i < _globalDegree.size(); i++)
  {
    notVisited.push_back(std::make_pair(i, _globalDegree[i]));
  }

  // Vector notVisited helps in running BFS
  // even when there are dijoind graphs
  while (!notVisited.empty())
  {
    label minNodeIndex = 0;

    for (label i = 0; i < notVisited.size(); i++)
    {
      if (notVisited[i].second < notVisited[minNodeIndex].second)
      {
        minNodeIndex = i;
      }
    }

    Q.push(notVisited[minNodeIndex].first);
    notVisited.erase(notVisited.begin() +  // minNodeIndex)
                     findIndex(notVisited, notVisited[Q.front()].first));
    // Simple BFS
    while (!Q.empty())
    {
      std::vector<label> *toSortPtr = new std::vector<label>();
      std::vector<label> &toSort = *toSortPtr;
      label y = Q.front();

      for (label i = 0; i < _nnz; i++)
      {
        if (_row[i] == Q.front() && _row[i] != _col[i] &&
            findIndex(notVisited, _col[i]) != -1)
        {
          toSort.push_back(_col[i]);
          notVisited.erase(notVisited.begin() + findIndex(notVisited, _col[i]));
        }
      }
      sort(toSort.begin(), toSort.end(), compareDegree);

      for (label i = 0; i < toSort.size(); i++)
      {
        Q.push(toSort[i]);
      }

      R.push_back(Q.front());
      Q.pop();
    }
  }

  DELETE_OBJECT_POINTER(Qptr);
  DELETE_OBJECT_POINTER(notVisitedPtr);
  POUT << "end CuthillMckee" << std::endl;
}

void ReorderingSSM::ReverseCuthillMckee()
{
  POUT << "in ReverseCuthillMckee" << std::endl;
  CuthillMckee();
  std::vector<label> &cuthill = *_vetexOrderPtr;

  label n = cuthill.size();

  if (n % 2 == 0)
  {
    n -= 1;
  }

  n = n / 2;

  //- reverse
  for (label i = 0; i <= n; i++)
  {
    label j = cuthill[cuthill.size() - 1 - i];
    cuthill[cuthill.size() - 1 - i] = cuthill[i];
    cuthill[i] = j;
  }

  //-
  label *ord = new label[_nRows];
  for (label i = 0; i < _nRows; i++)
  {
    ord[i] = i;
  }
  label otto;
  for (label i = 0; i < _nRows - 1; i++)
  {
    for (label j = _nRows - 1; j > i; j--)
    {
      if (cuthill[j - 1] > cuthill[j])
      {
        otto = cuthill[j - 1];
        cuthill[j - 1] = cuthill[j];
        cuthill[j] = otto;

        otto = ord[j - 1];
        ord[j - 1] = ord[j];
        ord[j] = otto;
      }
    }
  }
  for (label i = 0; i < _nRows; i++)
  {
    cuthill[i] = ord[i];
  }

  label Max = 0;
  for (label i = 0; i < _nnz; i++)
  {
    label newrow = cuthill[_row[i]];
    label newcol = cuthill[_col[i]];

    label otto = abs(newrow - newcol);
    if (otto > Max)
    {
      Max = otto;
    }
  }
  POUT << "Bandwidth after  reordering: " << Max << std::endl;
  DELETE_POINTER(ord);
  POUT << "end ReverseCuthillMckee" << std::endl;
}

std::vector<label> *ReorderingSSM::getVetexOrder() { return _vetexOrderPtr; }

std::vector<label> *ReorderingSSM::getEdgeOrder()
{
  // need to change
  std::vector<label> &edgeorder = *_edgeOrderPtr;

  label otto;
  label *newrow = new label[_nnz];
  label *newcol = new label[_nnz];
  label *order = new label[_nnz];

  for (label i = 0; i < _nnz; i++)
  {
    edgeorder.push_back(i);
    order[i] = i;
    newrow[i] = (*_vetexOrderPtr)[_row[i]];
    newcol[i] = (*_vetexOrderPtr)[_col[i]];
  }

  for (label i = 0; i < _nnz - 1; i++)
  {
    for (label j = _nnz - 1; j > i; j--)
    {
      if (newcol[j - 1] > newcol[j])
      {
        otto = newcol[j - 1];
        newcol[j - 1] = newcol[j];
        newcol[j] = otto;

        otto = newrow[j - 1];
        newrow[j - 1] = newrow[j];
        newrow[j] = otto;

        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;
      }
    }
  }
  for (label i = 0; i < _nnz - 1; i++)
  {
    for (label j = _nnz - 1; j > i; j--)
    {
      if (newrow[j - 1] > newrow[j])
      {
        otto = newcol[j - 1];
        newcol[j - 1] = newcol[j];
        newcol[j] = otto;
        otto = newrow[j - 1];
        newrow[j - 1] = newrow[j];
        newrow[j] = otto;

        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;
      }
    }
  }
  for (label i = 0; i < _nnz - 1; i++)
  {
    for (label j = _nnz - 1; j > i; j--)
    {
      if (order[j - 1] > order[j])
      {
        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;

        otto = edgeorder[j - 1];
        edgeorder[j - 1] = edgeorder[j];
        edgeorder[j] = otto;
      }
    }
  }
  DELETE_POINTER(newrow);
  DELETE_POINTER(newcol);
  DELETE_POINTER(order);
  return _edgeOrderPtr;
}

void rcmCOO_nowrite(const label nnz,
                    const label nRows,
                    const label *row,
                    const label *col,
                    label *postVetexOrder,
                    label *postValOrder)
{
  ReorderingSSM m(nnz, nRows, row, col);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);
  std::vector<label> *EdgeorderPtr = m.getEdgeOrder();
  std::vector<label> &Edgeorder = *(EdgeorderPtr);

  for (label ij = 0; ij < Cellorder.size(); ij++)
  {
    postVetexOrder[ij] = Cellorder[ij];
  }
  for (label ik = 0; ik < Edgeorder.size(); ik++)
  {
    postValOrder[ik] = Edgeorder[ik];
  }
}

void rcmCOO_rewrite(
    const label nnz, const label nRows, label *row, label *col, scalar *val)
{
  ReorderingSSM m(nnz, nRows, row, col);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);
  std::vector<label> *EdgeorderPtr = m.getEdgeOrder();
  std::vector<label> &Edgeorder = *(EdgeorderPtr);

  label *ncol = new label[nnz];
  label *nrow = new label[nnz];
  scalar *newval = new scalar[nnz];
  for (label j = 0; j < nnz; j++)
  {
    nrow[Edgeorder[j]] = Cellorder[row[j]];
    ncol[Edgeorder[j]] = Cellorder[col[j]];
    newval[Edgeorder[j]] = val[j];
  }
  for (label j = 0; j < nnz; j++)
  {
    col[j] = ncol[j];
    row[j] = nrow[j];
    val[j] = newval[j];
  }

  DELETE_POINTER(ncol);
  DELETE_POINTER(nrow);
  DELETE_POINTER(newval);
}

void rcmCSR_nowrite(const label nnz,
                    const label nRows,
                    const label *rowsOffset,
                    const label *col,
                    label *postVetexOrder,
                    label *postValOrder,
                    label *newRowsOffset)
{
  label *row = new label[nnz];
  label j = 0;
  for (label i = 0; i < nnz; i++)
  {
    if (i >= rowsOffset[j + 1])
    {
      j++;
    }
    row[i] = j;
  }

  ReorderingSSM m(nnz, nRows, row, col);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);
  std::vector<label> *EdgeorderPtr = m.getEdgeOrder();
  std::vector<label> &Edgeorder = *(EdgeorderPtr);

  for (label ij = 0; ij < Cellorder.size(); ij++)
  {
    postVetexOrder[ij] = Cellorder[ij];
  }
  for (label ik = 0; ik < Edgeorder.size(); ik++)
  {
    postValOrder[ik] = Edgeorder[ik];
  }
  newRowsOffset[0] = 0;

  for (label j = 1; j < nRows + 1; j++)
  {
    newRowsOffset[Cellorder[j - 1] + 1] = rowsOffset[j] - rowsOffset[j - 1];
  }
  for (label j = 1; j < nRows + 1; j++)
  {
    newRowsOffset[j] += newRowsOffset[j - 1];
  }

  DELETE_POINTER(row);
}

void rcmCSR_rewrite(const label nnz,
                    const label nRows,
                    label *rowsOffset,
                    label *col,
                    scalar *value)
{
  label *row = new label[nnz];
  label k = 0;
  for (label i = 0; i < nnz; i++)
  {
    if (i >= rowsOffset[k + 1])
    {
      k++;
    }
    row[i] = k;
  }
  ReorderingSSM m(nnz, nRows, row, col);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);
  std::vector<label> *EdgeorderPtr = m.getEdgeOrder();
  std::vector<label> &Edgeorder = *(EdgeorderPtr);

  label *newcol = new label[nnz];
  scalar *newval = new scalar[nnz];

  for (label j = 0; j < nnz; j++)
  {
    newcol[Edgeorder[j]] = Cellorder[col[j]];
    newval[Edgeorder[j]] = value[j];
  }

  for (label j = 0; j < nnz; j++)
  {
    col[j] = newcol[j];
    value[j] = newval[j];
  }

  label *newRowsOffset = new label[nRows + 1];
  newRowsOffset[0] = 0;
  for (label j = 1; j < nRows + 1; j++)
  {
    newRowsOffset[Cellorder[j - 1] + 1] = rowsOffset[j] - rowsOffset[j - 1];
  }
  for (label j = 1; j < nRows + 1; j++)
  {
    rowsOffset[j] = newRowsOffset[j] + newRowsOffset[j - 1];
  }

  DELETE_POINTER(row);
  DELETE_POINTER(newRowsOffset);
  DELETE_POINTER(newcol);
  DELETE_POINTER(newval);
}

void rcmLDU_nowrite(const label nnz,
                    const label nRows,
                    const label *row,
                    const label *col,
                    label *postVetexOrder,
                    label *postEdgeOrder)
{
  label len = 2 * nnz;
  label *col2 = new label[len];
  label *row2 = new label[len];
  label *nColsInRow = new label[nRows];
  label *colsOffset = new label[nRows + 1];
  colsOffset[0] = 0;
  for (label i = 0; i < nRows; i++)
  {
    nColsInRow[i] = 0;
    colsOffset[i + 1] = 0;
  }
  for (label i = 0; i < nnz; i++)
  {
    label r = row[i];
    label c = col[i];

    nColsInRow[r]++;
    nColsInRow[c]++;
  }

  for (label i = 0; i < nRows; i++)
  {
    colsOffset[i + 1] = colsOffset[i] + nColsInRow[i];
    nColsInRow[i] = 0;
  }

  for (label i = 0; i < nnz; i++)
  {
    label r = row[i];
    label c = col[i];

    row2[colsOffset[r] + nColsInRow[r]] = r;
    col2[colsOffset[r] + nColsInRow[r]] = c;
    nColsInRow[r]++;

    row2[colsOffset[c] + nColsInRow[c]] = c;
    col2[colsOffset[c] + nColsInRow[c]] = r;
    nColsInRow[c]++;
  }
  DELETE_POINTER(nColsInRow);
  DELETE_POINTER(colsOffset);

  ReorderingSSM m(len, nRows, row2, col2);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);
  label *order = new label[nnz];

  for (label it = 0; it < nRows; it++)
  {
    postVetexOrder[it] = Cellorder[it];
  }

  for (label it = 0; it < nnz; it++)
  {
    postEdgeOrder[it] = it;
    order[it] = it;
  }

  // need to change
  label *nrow = new label[nnz];
  label *ncol = new label[nnz];

  label otto;
  for (label i = 0; i < nnz; i++)
  {
    nrow[i] = Cellorder[row[i]];
    ncol[i] = Cellorder[col[i]];
    if (nrow[i] > ncol[i])
    {
      order[i] = -order[i] - 1;
      otto = nrow[i];
      nrow[i] = ncol[i];
      ncol[i] = otto;
    }
  }

  for (label i = 0; i < nnz; i++)
  {
  }

  for (label i = 0; i < nnz - 1; i++)
  {
    for (label j = nnz - 1; j > i; j--)
    {
      if (ncol[j - 1] > ncol[j])
      {
        otto = ncol[j - 1];
        ncol[j - 1] = ncol[j];
        ncol[j] = otto;

        otto = nrow[j - 1];
        nrow[j - 1] = nrow[j];
        nrow[j] = otto;

        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;
      }
    }
  }
  for (label i = 0; i < nnz - 1; i++)
  {
    for (label j = nnz - 1; j > i; j--)
    {
      if (nrow[j - 1] > nrow[j])
      {
        otto = ncol[j - 1];
        ncol[j - 1] = ncol[j];
        ncol[j] = otto;

        otto = nrow[j - 1];
        nrow[j - 1] = nrow[j];
        nrow[j] = otto;

        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;
      }
    }
  }

  for (label i = 0; i < nnz; i++)
  {
    if (order[i] < 0)
    {
      order[i] = -order[i] - 1;
      postEdgeOrder[i] = -postEdgeOrder[i] - 1;
    }
  }
  for (label i = 0; i < nnz - 1; i++)
  {
    for (label j = nnz - 1; j > i; j--)
    {
      if (order[j - 1] > order[j])
      {
        otto = order[j - 1];
        order[j - 1] = order[j];
        order[j] = otto;

        otto = postEdgeOrder[j - 1];
        postEdgeOrder[j - 1] = postEdgeOrder[j];
        postEdgeOrder[j] = otto;
      }
    }
  }

  DELETE_POINTER(ncol);
  DELETE_POINTER(nrow);
  DELETE_POINTER(col2);
  DELETE_POINTER(row2);
}

void rcmLDU_rewrite(const label nnz,
                    const label nRows,
                    label *row,
                    label *col,
                    scalar *diagVal,
                    scalar *upperVal,
                    scalar *lowerVal)
{
  label len = 2 * nnz;
  label *col2 = new label[len];
  label *row2 = new label[len];
  label *nColsInRow = new label[nRows];
  label *colsOffset = new label[nRows + 1];
  colsOffset[0] = 0;
  for (label i = 0; i < nRows; i++)
  {
    nColsInRow[i] = 0;
    colsOffset[i + 1] = 0;
  }
  for (label i = 0; i < nnz; i++)
  {
    label r = row[i];
    label c = col[i];

    nColsInRow[r]++;
    nColsInRow[c]++;
  }

  for (label i = 0; i < nRows; i++)
  {
    colsOffset[i + 1] = colsOffset[i] + nColsInRow[i];
    nColsInRow[i] = 0;
  }

  for (label i = 0; i < nnz; i++)
  {
    label r = row[i];
    label c = col[i];

    row2[colsOffset[r] + nColsInRow[r]] = r;
    col2[colsOffset[r] + nColsInRow[r]] = c;
    nColsInRow[r]++;

    row2[colsOffset[c] + nColsInRow[c]] = c;
    col2[colsOffset[c] + nColsInRow[c]] = r;
    nColsInRow[c]++;
  }
  DELETE_POINTER(nColsInRow);
  DELETE_POINTER(colsOffset);

  ReorderingSSM m(len, nRows, row2, col2);

  std::vector<label> *CellorderPtr = m.getVetexOrder();
  std::vector<label> &Cellorder = *(CellorderPtr);

  label temp;
  scalar tem;
  label *newcol = new label[nnz];
  label *newrow = new label[nnz];
  scalar *newdiag = new scalar[nRows];

  for (label i = 0; i < nRows; i++)
  {
    newdiag[Cellorder[i]] = diagVal[i];
  }
  for (label i = 0; i < nRows; i++)
  {
    diagVal[i] = newdiag[i];
  }
  for (label i = 0; i < nnz; i++)
  {
    newrow[i] = Cellorder[row[i]];
    newcol[i] = Cellorder[col[i]];
  }

  for (label j = 0; j < nnz; j++)
  {
    if (newrow[j] > newcol[j])
    {
      tem = upperVal[j];
      upperVal[j] = lowerVal[j];
      lowerVal[j] = tem;
      temp = newrow[j];
      newrow[j] = newcol[j];
      newcol[j] = temp;
    }
  }
  for (label j = 0; j < nnz; j++)
  {
    col[j] = newcol[j];
    row[j] = newrow[j];
  }

  label otto;
  scalar ott;
  for (label i = 0; i < nnz - 1; i++)
  {
    for (label j = nnz - 1; j > i; j--)
    {
      if (col[j - 1] > col[j])
      {
        otto = col[j - 1];
        col[j - 1] = col[j];
        col[j] = otto;

        otto = row[j - 1];
        row[j - 1] = row[j];
        row[j] = otto;

        ott = upperVal[j - 1];
        upperVal[j - 1] = upperVal[j];
        upperVal[j] = ott;
        ott = lowerVal[j - 1];
        lowerVal[j - 1] = lowerVal[j];
        lowerVal[j] = ott;
      }
    }
  }
  for (label i = 0; i < nnz - 1; i++)
  {
    for (label j = nnz - 1; j > i; j--)
    {
      if (row[j - 1] > row[j])
      {
        otto = col[j - 1];
        col[j - 1] = col[j];
        col[j] = otto;
        otto = row[j - 1];
        row[j - 1] = row[j];
        row[j] = otto;
        ott = upperVal[j - 1];
        upperVal[j - 1] = upperVal[j];
        upperVal[j] = ott;
        ott = lowerVal[j - 1];
        lowerVal[j - 1] = lowerVal[j];
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
