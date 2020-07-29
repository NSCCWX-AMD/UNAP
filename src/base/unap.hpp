/// \file unap.hpp
/// \brief brief information to be added
///
/// Detailed information to be added
///
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30

#ifndef UNAP_HPP
#define UNAP_HPP

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <typeinfo>

#include "utilities.h"
// #include "swMacro.h"

// define sw data type
typedef int swInt;
typedef int swInt32;
typedef long swInt64;

typedef double swFloat;
typedef float swFloat32;
typedef double swFloat64;

namespace UNAP
{
//- sleep for the specified number of seconds
unsigned int sleep(const unsigned int);

//- small scalar for the use in solvers
#define SMALL (1.0e-37)
#define VSMALL (1.0e-128)

//- large scalar for the use in solvers
#define GREAT (1.0e+37)
#define VGREAT (1.0e+128)

//- check if pointer exists
#define CHECK_POINTER(ptr)                                                 \
  if (!ptr)                                                                \
  {                                                                        \
    UNAPCOUT << "ERROR in " << __FILE__ << " " << __LINE__ << ": " << #ptr \
             << " is NULL!" << ENDL;                                       \
    ERROR_EXIT;                                                            \
  }

#define ALLOCATE_POINTER(ptr, oldObj, T) \
  {                                      \
    DELETE_OBJECT_POINTER(ptr)           \
    ptr = new T(oldObj.size());          \
    T &newObj = *ptr;                    \
    newObj = oldObj;                     \
  }

template <typename T>
T mag(T v)
{
  T absV;
  if (typeid(T) == typeid(label))
  {
    absV = abs(v);
  }
  else if (typeid(T) == typeid(scalar))
  {
    absV = fabs(v);
  }
  return absV;
}

#ifdef SW_SLAVE
#define IFNOT_SWACC if (nCells < accUsingSize)

#define IFNOT_SWACC_SMALL if (nCells < 2500)

#define SpMVAccSize 10000
#else
#define IFNOT_SWACC

#define IFNOT_SWACC_SMALL
#endif

#define printMessage(ss)       \
  MPI_Barrier(MPI_COMM_WORLD); \
  UNAPCOUT << ss << ENDL;

}  // namespace UNAP

#endif
