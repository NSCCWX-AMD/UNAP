#ifndef INTERFACES_HPP
#define INTERFACES_HPP

#include "PtrList.hpp"
#include "patch.hpp"

namespace UNAP
{
class patch;

class interfaces
{
 private:
  PtrList<patch> &patches_;

  //- send buffer for initializing and updating matrix interfaces
  mutable scalar *sendBuffer_;

  //- receive buffer for initializing and updating matrix interfaces
  mutable scalar *recvBuffer_;

  //- mpi requests of sends and receives for updating interfaces
  mutable MPI_Request *sendRecvRequests_;

  mutable string *sendTaskName_;

  mutable string *recvTaskName_;

  //- locPosition
  mutable label *locPosition_;

  //- destRank
  mutable label *destRank_;

  //- communicator
  mutable Communicator *commcator_;

 public:
  interfaces(const label size, Communicator *other_comm);

  interfaces(PtrList<patch> &patches, Communicator *other_comm);

  virtual ~interfaces();

  //- Initialize the update of interfaced interfaces
  //  for matrix operations
  virtual void initMatrixInterfaces(const scalarVector &psi) const;

  //- update interfaced interfaces for matrix operations
  virtual void updateMatrixInterfaces(scalarVector &result) const;

  virtual patch &patchList(const label i) const { return patches_[i]; }

  PtrList<patch> &patchList() { return patches_; }

  virtual label size() const { return patches_.size(); }

  void reorderIntFaceCells(const label *cellMap);

  // - set communicator
  void setCommunicator(Communicator *other_comm) { commcator_ = other_comm; }

  // - get communicator
  Communicator *getCommunicator() const { return commcator_; }
};

}  // namespace UNAP
#endif  //- INTERFACES_HPP
