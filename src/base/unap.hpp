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

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <typeinfo>
#include <math.h>

// #include "swMacro.h"

namespace UNAP
{
	#define COUT std::cout
	#define ENDL std::endl
	#define EXIT exit(0)
	#define ERROR_EXIT exit(1)
}


//- label
#if defined(LABEL_INT)
//- define label as a 32-bit int
namespace UNAP
{
	typedef int label;
}
#elif defined(LABEL_LONG)
//- define label as a 64-bit int
namespace UNAP
{
	typedef long int label;
}
#else
{
	COUT << "Error: label size is not defined." << ENDL;
	ERROR_EXIT;
}
#endif


//- scalar
#if defined(SCALAR_FLOAT)
// define scalar as a float
namespace UNAP
{
	typedef float scalar;
}
#elif defined(SCALAR_DOUBLE)
//- define scalar as a double
namespace UNAP
{
	typedef double scalar;
}
#else
{
	COUT << "Error: scalar size is not defined." << ENDL;
	ERROR_EXIT;
}
#endif

typedef int    swInt;
typedef int    swInt32;
typedef long   swInt64;

typedef double swFloat;
typedef float  swFloat32;
typedef double swFloat64;

namespace UNAP
{
	//- sleep for the specified number of seconds
	unsigned int sleep(const unsigned int);

	//- simplify loop
	#define forAll(i, length) for(int i=0; i<length; i++)

	//- small scalar for the use in solvers
	#define SMALL (1.0e-37)
    #define VSMALL (1.0e-128)

	//- large scalar for the use in solvers
	#define GREAT (1.0e+37)
    #define VGREAT (1.0e+128)

	//- check if pointer exists
	#define CHECK_POINTER(ptr) if(!ptr) \
    { \
    	COUT << "ERROR in " << __FILE__ << " " << __LINE__ \
			 << ": " << #ptr << " is NULL!" << ENDL; \
		ERROR_EXIT; \
    }

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

    #define ALLOCATE_POINTER(ptr, oldObj, T) \
    { \
        DELETE_OBJECT_POINTER(ptr) \
        ptr = new T(oldObj.size()); \
        T& newObj = *ptr; \
        newObj = oldObj; \
    }

    #define MAX(x, y) ( ((x) > (y)) ? (x) : (y) )
    #define MIN(x, y) ( ((x) < (y)) ? (x) : (y) )

    template<typename T>
    T mag(T v)
    {
        T absV;
        if(typeid(T) == typeid(label))
        {
            absV = abs(v);
        }
        else if(typeid(T) == typeid(scalar))
        {
            absV = fabs(v);
        }
        return absV;
    }

#ifdef SW_SLAVE
    #define IFNOT_SWACC \
    if(nCells < accUsingSize)

    #define IFNOT_SWACC_SMALL \
    if(nCells < 2500)

    #define SpMVAccSize 10000
#else
    #define IFNOT_SWACC

    #define IFNOT_SWACC_SMALL
#endif


#define printMessage(ss) \
MPI_Barrier(MPI_COMM_WORLD); \
if(!MYID) \
{ \
    COUT << ss << ENDL; \
}

}

#endif
