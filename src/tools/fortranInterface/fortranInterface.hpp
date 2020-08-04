#ifndef FORTRANINTERFACE_HPP
#define FORTRANINTERFACE_HPP

#include "unap.hpp"

namespace UNAP
{
#ifdef __cplusplus
extern "C"
{
#endif

  // 初始化MPI环境 返回communicator
  void comminit_(long int *CommPtr);

  void commgetmyidsize_(long int *CommPtr, int *rank, int *size);

  void ldumatrixcreat_(long int *APtrPtr,
                       label *nCells,
                       label *upperSize,
                       label *lowerAddr,
                       label *upperAddr,
                       scalar *lower,
                       scalar *diag,
                       scalar *upper,
                       long int *commPtr);

  void coo2ldumatrixcreat_(label64 *APtrPtr,
                           const scalar *dataPtr,
                           const label *rowsPtr,
                           const label *columnPtr,
                           const label *nCells,
                           const label *size,
                           const label *symm,
                           long int *commPtr);

  void csr2ldumatrixcreat_(label64 *APtrPtr,
                           const scalar *dataPtr,
                           const label *compRowsPtr,
                           const label *columnPtr,
                           const label *nCells,
                           const label *size,
                           const label *symm,
                           long int *commPtr);

  void matrixinterfacescreat_(label64 *APtrPtr,
                              const label *nNeiProcs,
                              const label *destRank,
                              const label *locPosition,
                              const label *faceCells,
                              const scalar *data);

  void pcgsolversolve_(scalar *x,
                       label64 *APtrPtr,
                       scalar *b,
                       label *nCells,
                       label *precond,
                       scalar *tol,
                       scalar *relTol,
                       label *maxIter,
                       label *minIter,
                       label *num_iterations,
                       scalar *final_res_norm);

  void pbicgstabsolversolve_(scalar *x,
                             label64 *APtrPtr,
                             scalar *b,
                             label *nCells,
                             label *precond,
                             scalar *tol,
                             scalar *relTol,
                             label *maxIter,
                             label *minIter,
                             label *num_iterations,
                             scalar *final_res_norm);

  void mgsolversolve_(scalar *x,
                      label64 *APtrPtr,
                      scalar *b,
                      label *nCells,
                      label *aggl,
                      label *smoother,
                      scalar *tol,
                      scalar *relTol,
                      label *maxIter,
                      label *minIter,
                      label *num_iterations,
                      scalar *final_res_norm,
                      scalar *faceAreaPtr);

  void reordercoo_(scalar *val,
                   label *row,
                   label *col,
                   label *nCellsPtr,
                   label *sizePtr,
                   long int *commPtr);

  void reorderuface__(label *row,
                      label *col,
                      label *nCellsPtr,
                      label *sizePtr,
                      label *newOrder,
                      long int *commPtr);

  void reorderlface__(label *row,
                      label *col,
                      label *nCellsPtr,
                      label *sizePtr,
                      label *newOrder,
                      long int *commPtr);

  void reordervalue__(scalar *val,
                      label *newOrder,
                      label *sizePtr,
                      long int *commPtr);

  void contruct_sw_matrix__(label64 *APtrPtr,
                            const label *nCellsPtr,
                            const label *rowsPtr,
                            const label *colsPtr,
                            const label *sizePtr,
                            long int *commPtr);

  void contruct_sw_matrix_interfaces__(label64 *APtrPtr,
                                       const label *nNeiProcsPtr,
                                       const label *destRankPtr,
                                       const label *offDiagRowsPtr,
                                       const label *offDiagStartsPtr);

#ifdef SW_SLAVE
  void construct_mlb_iterator__(label64 *APtrPtr);
#endif

  void fill_sw_matrix_coefficients__(label64 *APtrPtr,
                                     const scalar *diagPtr,
                                     const scalar *upperPtr,
                                     const scalar *lowerPtr);

  void fill_sw_matrix_interfaces_coefficients__(label64 *APtrPtr,
                                                const scalar *offDiagCoeffs);

  //- mg
  void contruct_solver_mg__(label64 *mgPtrPtr,
                            label64 *APtrPtr,
                            label64 *AgglPtr,
                            const scalar *faceAreaPtr,
                            const label *smTypePtr,
                            const label *maxLevelsPtr,
                            const label *nCellsCoarsestPtr);

#ifdef SW_SLAVE
  void mg_coarse_mlb__(label64 *agglPtr);
#endif

  void sw_solver_mg_set_maxiter__(label64 *solverPtrPtr,
                                  const label *maxIterPtr);

  void sw_solver_mg_set_miniter__(label64 *solverPtrPtr,
                                  const label *minIterPtr);

  void sw_solver_mg_set_tol__(label64 *solverPtrPtr, const scalar *tolPtr);

  void sw_solver_mg_set_reltol__(label64 *solverPtrPtr,
                                 const scalar *reltolPtr);

  void sw_solver_mg_set_npresweeps__(label64 *solverPtrPtr,
                                     const label *numPtr);

  void sw_solver_mg_set_npostsweeps__(label64 *solverPtrPtr,
                                      const label *numPtr);

  void sw_solver_mg_set_nfinestsweeps__(label64 *solverPtrPtr,
                                        const label *numPtr);

  void sw_solve_mg__(label64 *mgPtrPtr,
                     label64 *APtrPtr,
                     scalar *xPtr,
                     scalar *bPtr,
                     scalar *res_normPtr);

  //- pbicgstab

  void contruct_solver_pbicgstab__(long int *solverPtrPtr, long int *commPtr);

  void sw_solver_pbicgstab_set_maxiter__(label64 *solverPtrPtr,
                                         const label *maxIterPtr);

  void sw_solver_pbicgstab_set_miniter__(label64 *solverPtrPtr,
                                         const label *minIterPtr);

  void sw_solver_pbicgstab_set_tol__(label64 *solverPtrPtr,
                                     const scalar *tolPtr);

  void sw_solver_pbicgstab_set_reltol__(label64 *solverPtrPtr,
                                        const scalar *reltolPtr);

  void sw_solver_pbicgstab_set_precond__(label64 *solverPtrPtr,
                                         label64 *APtrPtr,
                                         const label *precondTypePtr);

  void sw_solve_pbicgstab__(label64 *solverPtrPtr,
                            label64 *APtrPtr,
                            scalar *xPtr,
                            scalar *bPtr,
                            scalar *res_normPtr);

  void sw_matrix_destroy__(label64 *APtrPtr);

  void sw_solver_destroy_mg__(label64 *solverPtrPtr);

  void sw_solver_destroy_pbicgstab__(label64 *solverPtrPtr);

  void createinterfaces__(label64 *APtrPtr,
                          label64 *offDiagRows,
                          label64 *offDiagCols,
                          label *offDiagPids,
                          label *cellNumPtr,
                          label *faceNumPtr,
                          label *postOrders);

  void printvector__(scalar *data, label *size, char *name, long int *commPtr);

#ifdef __cplusplus
}
#endif

}  // namespace UNAP

#endif  //- FORTRANINTERFACE_HPP
