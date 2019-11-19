#ifndef REORDERMATRIX_HPP
#define REORDERMATRIX_HPP

#include "lduMatrix.hpp"

namespace UNAP
{

void reorderLDUMatrix(lduMatrix& lduA);

void reorderLDUFaces(lduMatrix& lduA);

void reorderLDUCells()




} //- namespace UNAP
#endif //- REORDERMATRIX_HPP
