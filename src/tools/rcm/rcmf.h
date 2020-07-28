#ifndef RCMF_H_
#define RCMF_H_

#include "unap.hpp"

#ifdef __cplusplus
extern "C" {
#endif

void rcmCOO_nowrite(const label nnz,
					const label nRows,
					const label *row,
					const label *col,
					label *postVetexOrder,
					label *postValOrder);

void rcmCOO_rewrite(const label nnz,
					const label nRows,
					label *row,
					label *col,
					scalar *val);

void rcmCSR_nowrite(const label nnz,
					const label nRows,
					const label* rowsOffset,
					const label* col,
					label* postVetexOrder,
					label* postValOrder,
					label* newRowsOffset);

void rcmCSR_rewrite(const label nnz,
					const label nRows,
					label* rowsOffset,
					label* col,
					scalar* value);

void rcmLDU_nowrite(const label nnz,
					const label nRows,
					const label* row,
					const label* col,
					label* postVetexOrder,
					label* postEdgeOrder);

void rcmLDU_rewrite(const label nnz,
					const label nRows,
					label* row,
					label* col,
					scalar* diagVal,
					scalar* upperVal,
					scalar* lowerVal);

#ifdef __cplusplus
}
#endif

#endif
