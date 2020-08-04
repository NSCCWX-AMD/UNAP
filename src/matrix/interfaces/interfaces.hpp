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

public:
  interfaces(const label size);

  interfaces(PtrList<patch> &patches);

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
};

}  // namespace UNAP
#endif  //- INTERFACES_HPP
