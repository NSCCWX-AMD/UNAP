#ifndef PATCH_HPP
#define PATCH_HPP

#include "unapVector.hpp"

namespace UNAP
{

class patch
{
private:
	//- number of faces in this patch
	label size_;

	//- current processor number
	label myProcNo_;

	//- neighbor processor number
	label neighbProcNo_;

	//- face-cell addressing
	mutable labelField* faceCellsPtr_;

	//- coefficients
	mutable scalarField* patchCoeffsPtr_;

	//- face restriction addressing if needed
	mutable labelField* faceRestrictAddressingPtr_;

public:

	patch
	(
		label size,
		label myProcNo,
		label neighbProcNo
	);

	virtual ~patch();

	//- return neighbor processor number
	inline label neighbProcNo() const
	{
		return neighbProcNo_;
	}

	inline void neighbProcNo(const label i)
	{
		neighbProcNo_ = i;
	}

	//- return current processor number
	inline label myProcNo() const
	{
		return myProcNo_;
	}

	inline void myProcNo(const label i)
	{
		myProcNo_ = i;
	}

	//- return faceCells
	inline labelField& faceCells() const
	{
		CHECK_POINTER(faceCellsPtr_)
		return *faceCellsPtr_;
	}

	inline void faceCells(labelField& a) const
	{
		DELETE_OBJECT_POINTER(faceCellsPtr_)
		faceCellsPtr_ = &a;
	}

	//- return patch coefficients
	inline scalar patchCoeffs(const label faceI) const
	{
		return (*patchCoeffsPtr_)[faceI];
	}

	inline scalarField& patchCoeffs() const
	{
		CHECK_POINTER(patchCoeffsPtr_)
		return *patchCoeffsPtr_;
	}

	inline void patchCoeffs(scalarField& a)
	{
		DELETE_OBJECT_POINTER(patchCoeffsPtr_)
		patchCoeffsPtr_ = &a;
	}

	inline label size() const
	{
		return size_;
	}

	inline void size(const label i)
	{
		size_ = i;
	}

	inline void faceRestrictAddressing(labelField& a) const
	{
		DELETE_OBJECT_POINTER(faceRestrictAddressingPtr_)
		faceRestrictAddressingPtr_ = &a;
	}

	inline labelField& faceRestrictAddressing() const
	{
		CHECK_POINTER(faceRestrictAddressingPtr_)
		return *faceRestrictAddressingPtr_;
	}

	void reorderPatchFaceCells(const label* cellMap);
};

} //- end namespace UNAP
#endif //- PATCH_HPP
