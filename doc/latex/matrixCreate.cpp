//- LDU matrix create

label rowsSize = ?;  ///< rank of the matrix
label upperSize = ?; ///< number of non-zeros in the upper

label* rows = ?; ///< upper row index, length is upperSize
label* cols = ?; ///< upper column index, length is upperSize

scalar* upperVal = ?; ///< upper non-zeros values
scalar* lowerVal = ?; ///< lower non-zeros values
scalar* diagVal  = ?; ///< diagonal non-zeros values

labelField rowMap(rows, upperSize);
labelField colMap(cols, upperSize);

scalarField upper(upperVal, upperSize);
scalarField diag (diagVal, rowsSize);

/// if matrix is value-symmetric, lowerVal == upperVal
scalarField& lower = upper;
scalarField* lowerPtr = NULL;
///- else create space for lower
if(lowerVal != upperVal)
{
    lowerPtr = new scalarField(lowerVal, upperSize);
    lower = *lowerPtr;
}

/// create LDU matrix
lduMatrix* lduAPtr = new lduMatrix
(
    rowsSize,
    rowMap,
    colMap,
    lower,
    diag,
    upper
);

DELETE_OBJECT_POINTER(lowerPtr);

lduMatrix lduA = *lduAPtr;

///- create right hand side b and to-be-solved x
scalar* bVal = ?;
scalar* xVal = ?;

scalarField b (bVal, rowsSize);
scalarField x (xVal, rowsSize);
