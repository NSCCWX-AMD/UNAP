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

#include <mpi.h>

#include "unap.hpp"

namespace UNAP
{
class unapMPI
{
private:
  static int myProcNo_;
  static int nProcs_;
  static bool parRun_;
  static CommData unapLabel_;
  static CommData unapScalar_;
  static Communicator unapCommunicator_;

public:
  static void initMPI();
  static void initMPI(Communicator &otherCommunicator);

  static void exitMPI();

  //- is this a parallel run?
  static bool &parRun() { return parRun_; }

  //- number of processes in parallel run
  static label nProcs() { return nProcs_; }

  //- number of this process
  static label myProcNo() { return myProcNo_; }

  static Communicator &unapCommunicator() { return unapCommunicator_; }

  //- int type
  static MPI_Datatype &unapLabel() { return unapLabel_; }

  //- double type
  static MPI_Datatype &unapScalar() { return unapScalar_; }
};

#define MYID unapMPI::myProcNo()
#define PARRUN unapMPI::parRun()
#define UNAPMPI_LABEL unapMPI::unapLabel()
#define UNAPMPI_SCALAR unapMPI::unapScalar()
#define NPROCS unapMPI::nProcs()

#define UNAPCOUT \
  if (!MYID) std::cout

template <typename T>
void reduceSum(T *v, label n = 1)
{
  if (PARRUN)
  {
    T vLocal[n];
    forAll(i, n) { vLocal[i] = v[i]; }

    CommData myType;
    if (typeid(T) == typeid(label))
    {
      myType = UNAPMPI_LABEL;
    }
    else if (typeid(T) == typeid(scalar))
    {
      myType = UNAPMPI_SCALAR;
    }
    UNAP::unapMPI::unapCommunicator().allReduce(
        "sum", &vLocal[0], v, n, myType, COMM_SUM);
    UNAP::unapMPI::unapCommunicator().finishTask("sum");
  }
}

}  // namespace UNAP

#endif  //- UNAPMPI_HPP
