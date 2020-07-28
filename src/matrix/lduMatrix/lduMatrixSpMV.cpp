#include "lduMatrix.hpp"

#ifdef SW_SLAVE
#include "wrappedInterface.h"
#endif

#ifdef SWTIMER
#include "swTimer.hpp"
#endif

#ifdef SW_SLAVE
	extern "C" define_e2v_FunPtr(swSpMV);
	extern "C" define_e2v_slaveFunPtr(swSpMV);
#endif

void UNAP::lduMatrix::spMV
(
	scalarVector& Apsi,
	const scalarVector& psi
) const
{
	this->initInterfaces(psi);

	scalar* __restrict__ ApsiPtr = Apsi.begin();
	const scalar* psiPtr = psi.begin();

	const scalar* diagPtr = this->diag().begin();

	const label* uPtr = this->upperAddr().begin();
	const label* lPtr = this->lowerAddr().begin();

	const scalar* upperPtr = this->upper().begin();
	const scalar* lowerPtr = this->lower().begin();

	const label nCells = this->nCells();
	const label nFaces = this->upper().size();

	forAll(cellI, nCells)
	{
		ApsiPtr[cellI] = 0.0;
	}

#ifdef SWTIMER
    swTimer::startTimer("SpMV");
#endif

#ifdef SW_SLAVE
    if(unatIter_)
    {
    	Arrays backEdgeData, frontEdgeData, selfConnData, vertexData, paraData;
		Arrays backEdgeData_p, frontEdgeData_p, selfConnData_p, vertexData_p;
		constructSingleArray(backEdgeData, 1, nFaces, COPYIN,
					(swFloat*)lowerPtr);
		constructSingleArray(frontEdgeData, 1, nFaces, COPYIN,
					(swFloat*)upperPtr);
		constructSingleArray(selfConnData, 1, nCells, COPYIN,
					(swFloat*)diagPtr);
		constructSingleArray(vertexData, 1, nCells, COPYIN,
					(swFloat*)psiPtr);
		addSingleArray(vertexData, 1, nCells, COPYOUT,
					(swFloat*)ApsiPtr);
		constructEmptyArray(backEdgeData_p);
		constructEmptyArray(frontEdgeData_p);
		constructEmptyArray(selfConnData_p);
		constructEmptyArray(vertexData_p);
		constructEmptyArray(paraData);
		FieldData data = {&backEdgeData, &frontEdgeData, &selfConnData,
			&vertexData};
		FieldData data_p = {&backEdgeData_p, &frontEdgeData_p, &selfConnData_p,
			&vertexData_p};
		coupledOperator cOpt;
		cOpt.fun_host  = swSpMV;
		cOpt.fun_slave = slave_swSpMV;
		cOpt.data = &data;
		cOpt.data_p = &data_p;

		// Compute
		this->unatIter_->edge2VertexIteration(&paraData, &cOpt, 1);
    }
    else
#endif
    {
    	for(label cell=0; cell<nCells; cell++)
	    {
	        ApsiPtr[cell] += diagPtr[cell]*psiPtr[cell];
	    }
	    for(label face=0; face<nFaces; face++)
	    {
	        ApsiPtr[uPtr[face]] += lowerPtr[face]*psiPtr[lPtr[face]];
	        ApsiPtr[lPtr[face]] += upperPtr[face]*psiPtr[uPtr[face]];
	    }
    }

#ifdef SWTIMER
    swTimer::endTimer("SpMV");
#endif
    this->updateInterfaces(Apsi);
}
