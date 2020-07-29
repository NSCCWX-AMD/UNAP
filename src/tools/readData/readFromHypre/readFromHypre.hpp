#ifndef READFROMHYPRE_HPP
#define READFROMHYPRE_HPP

#include <map>

#include "lduMatrix.hpp"

namespace UNAP
{
void constructLDUMatrixFromHypre(lduMatrix &lduA, const char *fileName);

void sortInterFaces(scalarVector &val,
                    labelVector &row,
                    labelVector &col,
                    const labelVector &faceToProcNO,
                    const label faceSize,
                    const labelVector &faceStart,
                    std::map<int, int> &mapProcNO,
                    const labelVector &neiProcNo,
                    const label procSize,
                    const labelVector &globalRowStart,
                    const labelVector &globalRowEnd);

void constructLDUInterfacesFromHypre(lduMatrix &lduA,
                                     const label nNeiProcs,
                                     const labelVector &destRank,
                                     const labelVector &locPosition,
                                     const labelVector &faceCells,
                                     const scalarVector &data);

void constructVectorFromHypre(scalarVector &b, const char *fileName);
}  // namespace UNAP
#endif
