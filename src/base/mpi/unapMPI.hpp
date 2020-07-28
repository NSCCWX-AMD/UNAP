/// \file unapMPI.hpp
/// \brief brief information to be added
///
/// Detailed information to be added
///
/// \author Hanfeng GU
/// \version 1.0
/// \date 2019-01-30

#ifndef UNAPMPI_HPP
#define UNAPMPI_HPP

#include "unap.hpp"
#include <mpi.h>

namespace UNAP
{

class unapMPI
{
private:
	static int myProcNo_;
	static int nProcs_;
	static bool parRun_;
	static MPI_Datatype unapLabel_;
	static MPI_Datatype unapScalar_;

public:

	static void initMPI();

    static void exitMPI();

	//- is this a parallel run?
    static bool& parRun()
    {
        return parRun_;
    }

	//- number of processes in parallel run
    static label nProcs()
    {
        return nProcs_;
    }

    //- number of this process
    static label myProcNo()
    {
        return myProcNo_;
    }

    //- int type
    static MPI_Datatype& unapLabel()
    {
    	return unapLabel_;
    }

    //- double type
    static MPI_Datatype& unapScalar()
    {
    	return unapScalar_;
    }
};


#define MYID       unapMPI::myProcNo()
#define PARRUN     unapMPI::parRun()
#define UNAPMPI_LABEL  unapMPI::unapLabel()
#define UNAPMPI_SCALAR unapMPI::unapScalar()
#define NPROCS     unapMPI::nProcs()


template<typename T>
void reduceSum(T& v)
{
	if(PARRUN)
	{
		T vLocal = v;
		MPI_Datatype myType;
		if(typeid(T) == typeid(label))
		{
			myType = UNAPMPI_LABEL;
		}
		else if(typeid(T) == typeid(scalar))
		{
			myType = UNAPMPI_SCALAR;
		}

        // MPI_Barrier(MPI_COMM_WORLD);
		MPI_Allreduce(&vLocal, &v, 1, myType, MPI_SUM, MPI_COMM_WORLD);
	}
}


template<typename T>
void reduceSum(T* v, label n)
{
    if(PARRUN)
    {
        T vLocal[n];
        forAll(i, n)
        {
            vLocal[i] = v[i];
        }

        MPI_Datatype myType;
        if(typeid(T) == typeid(label))
        {
            myType = UNAPMPI_LABEL;
        }
        else if(typeid(T) == typeid(scalar))
        {
            myType = UNAPMPI_SCALAR;
        }

        // MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&vLocal[0], v, n, myType, MPI_SUM, MPI_COMM_WORLD);
    }
}

} //- end namespace UNAP

#endif //- UNAPMPI_HPP
