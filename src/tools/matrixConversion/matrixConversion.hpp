#ifndef MATRIXCONVERSION_HPP
#define MATRIXCONVERSION_HPP

#include "lduMatrix.hpp"

namespace UNAP
{
void sortData(scalarVector &data,
              const labelVector &order,
              const labelVector &cellFaces,
              Communicator *other_comm);

//- reorder COO matrix, to make it ordered in columns in every row
//- to use this, the input COO matrix must be ordered in rows
//- row, col must be started from 0
void reorderCOO(scalar *dataPtr,
                label *rowsPtr,
                label *columnPtr,
                const label nCells,
                const label size,
                Communicator *other_comm);

void reorderValue(scalar *val,
                  const label *newOrder,
                  const label size,
                  Communicator *other_comm);

void reorderUFace(label *rowsPtr,
                  label *columnPtr,
                  const label nCells,
                  const label size,
                  label *newOrder,
                  Communicator *other_comm);

void reorderLFace(label *rowsPtr,
                  label *columnPtr,
                  const label nCells,
                  const label size,
                  label *newOrder,
                  Communicator *other_comm);

//- only for diagonal part
lduMatrix &coo2ldu(const scalar *dataPtr,
                   const label *rowsPtr,
                   const label *columnPtr,
                   const label nCells,
                   const label size,
                   const bool symm,  //- symm refers to the data
                   Communicator *other_comm);

//- only for diagonal part
lduMatrix &csr2ldu(const scalar *dataPtr,
                   const label *compRowsPtr,
                   const label *columnPtr,
                   const label nCells,
                   const label size,
                   const bool symm,  //- symm refers to the data
                   Communicator *other_comm);

}  // namespace UNAP

#endif  //- MATRIXCONVERSION_HPP
