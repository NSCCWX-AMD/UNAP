#include "unapMPI.hpp"

CommData UNAP::unapMPI::unapLabel_ = COMM_INT;
CommData UNAP::unapMPI::unapScalar_ = COMM_DOUBLE;

UNAP::unapMPI::unapMPI() : unapCommunicator_(NULL)
{
  int initialized;
  MPI_Initialized(&initialized);
  if (!initialized)
  {
    COMM::init(NULL, NULL);
    unapCommunicator_ = &COMM::getGlobalComm();
  }

  nProcs_ = unapCommunicator_->getMySize();
  myProcNo_ = unapCommunicator_->getMyId();

  if (nProcs_ > 1)
  {
    parRun_ = true;
  }
  else
  {
    parRun_ = false;
  }

  if (typeid(label) == typeid(int))
  {
    unapLabel_ = COMM_INT;
  }
  else if (typeid(label) == typeid(long int))
  {
    unapLabel_ = COMM_LONG;
  }
  else
  {
    if (!myProcNo_)
      std::cout << "Error: label is not neither a int nor a long int type!"
                << ENDL;

    ERROR_EXIT;
  }

  if (typeid(scalar) == typeid(double))
  {
    unapScalar_ = COMM_DOUBLE;
  }
  else if (typeid(scalar) == typeid(float))
  {
    unapScalar_ = COMM_FLOAT;
  }
  else
  {
    if (!myProcNo_)
      std::cout << "Error: scalar is not neither a float nor a double type!"
                << ENDL;

    ERROR_EXIT;
  }
}

void UNAP::unapMPI::initMPI(Communicator *other_comm)
{
  unapCommunicator_ = other_comm;

  nProcs_ = unapCommunicator_->getMySize();
  myProcNo_ = unapCommunicator_->getMyId();

  if (nProcs_ > 1)
  {
    parRun_ = true;
  }
  else
  {
    parRun_ = false;
  }

  if (typeid(label) == typeid(int))
  {
    unapLabel_ = COMM_INT;
  }
  else if (typeid(label) == typeid(long int))
  {
    unapLabel_ = COMM_LONG;
  }
  else
  {
    if (!myProcNo_)
      std::cout << "Error: label is not neither a int nor a long int type!"
                << ENDL;
    ERROR_EXIT;
  }

  if (typeid(scalar) == typeid(double))
  {
    unapScalar_ = COMM_DOUBLE;
  }
  else if (typeid(scalar) == typeid(float))
  {
    unapScalar_ = COMM_FLOAT;
  }
  else
  {
    if (!myProcNo_)
      std::cout << "Error: scalar is not neither a float nor a double type!"
                << ENDL;
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
