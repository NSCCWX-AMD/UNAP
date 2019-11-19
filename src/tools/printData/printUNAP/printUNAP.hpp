#ifndef PRINTUNAP_HPP
#define PRINTUNAP_HPP

#include "lduMatrix.hpp"

#include <fstream>
#include <sstream>
#include <iomanip>

#define FILEOPEN(fout, fileName) \
std::ostringstream os; \
os << fileName << "_" << UNAP::MYID << ".txt"; \
fout.open(os.str().c_str())

#define FILECLOSE(fout, fileName) fout.close()

namespace UNAP
{
	void printLDUMatrix(const lduMatrix& A, const char* name);

	template<typename T>
	void printVector(const T& b, const char* name);

	void printInterfaces(const lduMatrix& A, const char* name);
}

template<typename T>
void UNAP::printVector(const T& v, const char* fileName)
{
	const label size = v.size();
    std::ofstream fout;

	FILEOPEN(fout, fileName);

    fout << "nCells: " << size << std::endl;
    forAll(i, size)
    {
    	fout << std::setiosflags(std::ios::scientific) << std::setprecision(15) << v[i] << std::endl;
    }

    FILECLOSE(fout, fileName);
}


#endif //- PRINTUNAP_HPP
