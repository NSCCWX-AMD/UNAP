#include <math.h>
#include "MG.hpp"
#include "lduAgglomeration.hpp"
#include "chebySmoother.hpp"
#include "lduGaussSeidelSmoother.hpp"
#include "matrixConversion.hpp"
#include "PBiCGStab.hpp"
#include "PCG.hpp"
#include "lduDiagPrecond.hpp"
#include "lduDICPrecond.hpp"
#include "lduDILUPrecond.hpp"
#include <sstream>
#include <string.h>
#include "printUNAP.hpp"
#include "readFromHypre.hpp"

#ifdef SW_SLAVE
#include "swAthread.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

//- test using data output from Hypre

using namespace UNAP;

#define LOCATEFILE(newName, fileName, dir) \
{ \
	std::ostringstream os; \
	os << dir << NPROCS << "/" << fileName << ".txt." << std::setfill('0') << std::setw(5) << MYID; \
	strcpy(newName, os.str().c_str()); \
} \


int main()
{
	/* Initialize MPI */
   	unapMPI::initMPI();
#ifdef SW_SLAVE
	swacc_init();
#endif

	const char* dir = "./exData/hypre/int_test/p";
   	char fileName[200];

   	if(PARRUN)
   	{
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

   	if(!MYID)

	{
		COUT << "Start reading data" << ENDL;
	}

   	lduMatrix lduA;
   	LOCATEFILE(fileName, "A_u", dir);
   	constructLDUMatrixFromHypre(lduA, fileName);

   	label nCells = lduA.size();
   	scalarField b(nCells);

   	LOCATEFILE(fileName, "b_u", dir);
   	constructVectorFromHypre(b, fileName);

   	if(PARRUN)
   	{
   		MPI_Barrier(MPI_COMM_WORLD);
   	}

   	if(!MYID)
	{
		COUT << "Finish reading data" << ENDL;
	}





#ifdef SW_SLAVE
	swacc_end();
#endif

#ifdef SWTIMER
	swTimer::printTimer();
#endif

	MPI_Finalize();
}
