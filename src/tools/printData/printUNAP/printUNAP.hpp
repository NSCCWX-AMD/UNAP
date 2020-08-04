#ifndef PRINTUNAP_HPP
#define PRINTUNAP_HPP

#include <fstream>
#include <iomanip>
#include <sstream>

#include "lduMatrix.hpp"

#define FILEOPEN(fout, fileName, MYID)     \
  std::ostringstream os;                   \
  os << fileName << "_" << MYID << ".txt"; \
  fout.open(os.str().c_str())

#define FILECLOSE(fout, fileName) fout.close()

namespace UNAP
{
void printLDUMatrix(const lduMatrix &A, const char *name);

template <typename T>
void printVector(const T &b, const char *name);

void printInterfaces(const lduMatrix &A, const char *name);
}  // namespace UNAP

template <typename T>
void UNAP::printVector(const T &v, const char *fileName)
{
  const label size = v.size();
  std::ofstream fout;

  FILEOPEN(fout, fileName, v.getCommunicator()->getMyId());

  fout << "nCells: " << size << std::endl;
  forAll(i, size)
  {
    fout << std::setiosflags(std::ios::scientific) << std::setprecision(15)
         << v[i] << std::endl;
  }

  FILECLOSE(fout, fileName);
}

#endif  //- PRINTUNAP_HPP
