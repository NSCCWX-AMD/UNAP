#include "unapMPI.hpp"

UNAP::label UNAP::unapMPI::myProcNo_ = 0;
UNAP::label UNAP::unapMPI::nProcs_ = 1;
bool  UNAP::unapMPI::parRun_ = false;
MPI_Datatype UNAP::unapMPI::unapLabel_ = MPI_INT;
MPI_Datatype UNAP::unapMPI::unapScalar_ = MPI_DOUBLE;

void UNAP::unapMPI::initMPI()
{
    int initialized;
    MPI_Initialized(&initialized);
    if(!initialized)
    {
        MPI_Init(NULL, NULL);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &nProcs_);
    MPI_Comm_rank(MPI_COMM_WORLD, &myProcNo_);

    if(nProcs_ > 1)
    {
    	parRun_ = true;
    }
    else
    {
    	parRun_ = false;
    }

    if(typeid(label) == typeid(int))
    {
    	unapLabel_ = MPI_INT;
    }
    else if(typeid(label) == typeid(long int))
    {
    	unapLabel_ = MPI_LONG;
    }
    else
    {
    	COUT << "Error: label is not neither a int nor a long int type!" << ENDL;
    	ERROR_EXIT;
    }

    if(typeid(scalar) == typeid(double))
    {
    	unapScalar_ = MPI_DOUBLE;
    }
    else if(typeid(scalar) == typeid(float))
    {
    	unapScalar_ = MPI_FLOAT;
    }
    else
    {
    	COUT << "Error: scalar is not neither a float nor a double type!" << ENDL;
    	ERROR_EXIT;
    }
}


void UNAP::unapMPI::exitMPI()
{
    int finalized;
    MPI_Finalized(&finalized);
    if (!finalized)
    {
        MPI_Finalize();
    }
}
