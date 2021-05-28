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
  int myProcNo_;
  int nProcs_;
  bool parRun_;
  static CommData unapLabel_;
  static CommData unapScalar_;
  Communicator *commcator_;

 public:
  unapMPI();

  void initMPI();
  void initMPI(Communicator *other_comm);

  void exitMPI();

  //- is this a parallel run?
  bool &parRun() { return parRun_; }

  //- number of processes in parallel run
  label nProcs() { return nProcs_; }

  //- number of this process
  label myProcNo() { return myProcNo_; }

  Communicator *unapCommunicator() { return commcator_; }

  //- int type
  static MPI_Datatype &unapLabel() { return unapLabel_; }

  //- double type
  static MPI_Datatype &unapScalar() { return unapScalar_; }
};

#define UNAPMPI_LABEL unapMPI::unapLabel()
#define UNAPMPI_SCALAR unapMPI::unapScalar()

template <typename T>
void reduceSum(T *v, Communicator *commcator, label n = 1)
{
  if (commcator->getMySize() > 1)
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
    commcator->allReduce("sum", &vLocal[0], v, n, myType, COMM_SUM);
    commcator->finishTask("sum");
  }
}

}  // namespace UNAP

#endif  //- UNAPMPI_HPP
