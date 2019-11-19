#ifndef READFROMOPENFOAM_HPP
#define READFROMOPENFOAM_HPP

#include "lduMatrix.hpp"

namespace UNAP
{

void constructLDUMatrixFromOpenFOAM(lduMatrix& lduA, const char* fileName);

void constructVectorFromOpenFOAM(scalarField& b, const char* fileName);

void constructLDUInterfacesFromOpenFOAM(lduMatrix& lduA, const char* fileName);

}
#endif
