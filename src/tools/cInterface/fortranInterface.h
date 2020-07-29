#ifndef FORTRANINTERFACE_H
#define FORTRANINTERFACE_H

#ifdef __cplusplus
extern "C"
{
#endif

  void ldumatrixcreat_(long int *APtrPtr,
                       int *nCells,
                       int *upperSize,
                       int *lowerAddr,
                       int *upperAddr,
                       double *lower,
                       double *diag,
                       double *upper);

  void coo2ldumatrixcreat_(long int *APtrPtr,
                           const double *dataPtr,
                           const int *rowsPtr,
                           const int *columnPtr,
                           const int *nCells,
                           const int *size,
                           const int *symm);

  void csr2ldumatrixcreat_(long int *APtrPtr,
                           const double *dataPtr,
                           const int *compRowsPtr,
                           const int *columnPtr,
                           const int *nCells,
                           const int *size,
                           const int *symm);

  void matrixinterfacescreat_(long int *APtrPtr,
                              const int *nNeiProcs,
                              const int *destRank,
                              const int *locPosition,
                              const int *faceCells,
                              const double *data);

  void pcgsolversolve_(double *x,
                       long int *APtrPtr,
                       double *b,
                       int *nCells,
                       int *precond,
                       double *tol,
                       double *relTol,
                       int *maxIter,
                       int *minIter,
                       int *num_iterations,
                       double *final_res_norm);

  void pbicgstabsolversolve_(double *x,
                             long int *APtrPtr,
                             double *b,
                             int *nCells,
                             int *precond,
                             double *tol,
                             double *relTol,
                             int *maxIter,
                             int *minIter,
                             int *num_iterations,
                             double *final_res_norm);

  void mgsolversolve_(double *x,
                      long int *APtrPtr,
                      double *b,
                      int *nCells,
                      int *aggl,
                      int *smoother,
                      double *tol,
                      double *relTol,
                      int *maxIter,
                      int *minIter,
                      int *num_iterations,
                      double *final_res_norm,
                      double *faceAreaPtr);

#ifdef __cplusplus
}
#endif

#endif  //- FORTRANINTERFACE_H
