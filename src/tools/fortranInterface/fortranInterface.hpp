#ifndef FORTRANINTERFACE_HPP
#define FORTRANINTERFACE_HPP

#include "unap.hpp"

namespace UNAP
{
#ifdef __cplusplus
extern "C"
{
#endif

  void ldumatrixcreat_(long int *APtrPtr,
                       label *nCells,
                       label *upperSize,
                       label *lowerAddr,
                       label *upperAddr,
                       scalar *lower,
                       scalar *diag,
                       scalar *upper);

  void coo2ldumatrixcreat_(long int *APtrPtr,
                           const scalar *dataPtr,
                           const label *rowsPtr,
                           const label *columnPtr,
                           const label *nCells,
                           const label *size,
                           const label *symm);

  void csr2ldumatrixcreat_(long int *APtrPtr,
                           const scalar *dataPtr,
                           const label *compRowsPtr,
                           const label *columnPtr,
                           const label *nCells,
                           const label *size,
                           const label *symm);

  void matrixinterfacescreat_(long int *APtrPtr,
                              const label *nNeiProcs,
                              const label *destRank,
                              const label *locPosition,
                              const label *faceCells,
                              const scalar *data);

  void pcgsolversolve_(scalar *x,
                       long int *APtrPtr,
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
                             long int *APtrPtr,
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
                      long int *APtrPtr,
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

  void reordercoo_(
      scalar *val, label *row, label *col, label *nCellsPtr, label *sizePtr);

  void reorderuface__(label *row,
                      label *col,
                      label *nCellsPtr,
                      label *sizePtr,
                      label *newOrder);

  void reorderlface__(label *row,
                      label *col,
                      label *nCellsPtr,
                      label *sizePtr,
                      label *newOrder);

  void reordervalue__(scalar *val, label *newOrder, label *sizePtr);

  void contruct_sw_matrix__(long int *APtrPtr,
                            const label *nCellsPtr,
                            const label *rowsPtr,
                            const label *colsPtr,
                            const label *sizePtr);

  void contruct_sw_matrix_interfaces__(long int *APtrPtr,
                                       const label *nNeiProcsPtr,
                                       const label *destRankPtr,
                                       const label *offDiagRowsPtr,
                                       const label *offDiagStartsPtr);

#ifdef SW_SLAVE
  void construct_mlb_iterator__(long int *APtrPtr);
#endif

  void fill_sw_matrix_coefficients__(long int *APtrPtr,
                                     const scalar *diagPtr,
                                     const scalar *upperPtr,
                                     const scalar *lowerPtr);

  void fill_sw_matrix_interfaces_coefficients__(long int *APtrPtr,
                                                const label *offDiagStartsPtr,
                                                const scalar *offDiagCoeffs);

  //- mg
  void contruct_solver_mg__(long int *mgPtrPtr,
                            long int *APtrPtr,
                            long int *AgglPtr,
                            const scalar *faceAreaPtr,
                            const label *smTypePtr,
                            const label *maxLevelsPtr,
                            const label *nCellsCoarsestPtr);

#ifdef SW_SLAVE
  void mg_coarse_mlb__(long int *agglPtr);
#endif

  void sw_solver_mg_set_maxiter__(long int *solverPtrPtr,
                                  const label *maxIterPtr);

  void sw_solver_mg_set_miniter__(long int *solverPtrPtr,
                                  const label *minIterPtr);

  void sw_solver_mg_set_tol__(long int *solverPtrPtr, const scalar *tolPtr);

  void sw_solver_mg_set_reltol__(long int *solverPtrPtr,
                                 const scalar *reltolPtr);

  void sw_solver_mg_set_npresweeps__(long int *solverPtrPtr,
                                     const label *numPtr);

  void sw_solver_mg_set_npostsweeps__(long int *solverPtrPtr,
                                      const label *numPtr);

  void sw_solver_mg_set_nfinestsweeps__(long int *solverPtrPtr,
                                        const label *numPtr);

  void sw_solve_mg__(long int *mgPtrPtr,
                     long int *APtrPtr,
                     scalar *xPtr,
                     scalar *bPtr,
                     scalar *res_normPtr);

  //- pbicgstab

  void contruct_solver_pbicgstab__(long int *solverPtrPtr);

  void sw_solver_pbicgstab_set_maxiter__(long int *solverPtrPtr,
                                         const label *maxIterPtr);

  void sw_solver_pbicgstab_set_miniter__(long int *solverPtrPtr,
                                         const label *minIterPtr);

  void sw_solver_pbicgstab_set_tol__(long int *solverPtrPtr,
                                     const scalar *tolPtr);

  void sw_solver_pbicgstab_set_reltol__(long int *solverPtrPtr,
                                        const scalar *reltolPtr);

  void sw_solver_pbicgstab_set_precond__(long int *solverPtrPtr,
                                         long int *APtrPtr,
                                         const label *precondTypePtr);

  void sw_solve_pbicgstab__(long int *solverPtrPtr,
                            long int *APtrPtr,
                            scalar *xPtr,
                            scalar *bPtr,
                            scalar *res_normPtr);

  void sw_matrix_destroy__(long int *APtrPtr);

  void sw_solver_destroy_mg__(long int *solverPtrPtr);

  void sw_solver_destroy_pbicgstab__(long int *solverPtrPtr);

  void createinterfaces__(label64 *APtrPtr,
                          label64 *offDiagRows,
                          label64 *offDiagCols,
                          label *offDiagPids,
                          label *cellNumPtr,
                          label *faceNumPtr,
                          label *postOrders);

  void printvector__(scalar *data, label *size, char *name);

#ifdef __cplusplus
}
#endif

}  // namespace UNAP

#endif  //- FORTRANINTERFACE_HPP
