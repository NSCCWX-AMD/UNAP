#ifndef RCM_H_
#define RCM_H_

#ifdef __cplusplus
extern "C" {
#endif

void rcmCOO_nowrite(const int nnz,
					const int nRows,
					const int *row,
					const int *col,
					int *postVetexOrder,
					int *postValOrder);

void rcmCOO_rewrite(const int nnz,
					const int nRows,
					int *row,
					int *col,
					double *val);

void rcmCSR_nowrite(const int nnz,
					const int nRows,
					const int* rowsOffset,
					const int* col,
					int* postVetexOrder,
					int* postValOrder,
					int* newRowsOffset);

void rcmCSR_rewrite(const int nnz,
					const int nRows,
					int* rowsOffset,
					int* col,
					double* value);

void rcmLDU_nowrite(const int nnz,
					const int nRows,
					const int* row,
					const int* col,
					int* postVetexOrder,
					int* postEdgeOrder);

void rcmLDU_rewrite(const int nnz,
					const int nRows,
					int* row,
					int* col,
					double* diagVal,
					double* upperVal,
					double* lowerVal);

#ifdef __cplusplus
}
#endif

#endif
