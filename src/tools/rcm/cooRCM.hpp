// C++ program for Implementation of
// Reverse Cuthill Mckee Algorithm
#ifndef cooRCM_HPP_
#define cooRCM_HPP_

#include <vector>
#include "unap.hpp"


namespace COORCM
{
class ReorderingSSM {
private:
	const label _nnz;    // number of non-zeros
	const label _nRows;  // number of rows
	const label* _row;   // array of rows
	const label* _col;   // array of columns

	static std::vector<label> _globalDegree;
	std::vector<label>* _vetexOrderPtr;
	std::vector<label>* _edgeOrderPtr;

	// label* _rowsOffset;

	void degreeGenerator();
    static bool compareDegree(label i, label j);

	void CuthillMckee();

	// Implementation of reverse Cuthill-Mckee algorithm
	void ReverseCuthillMckee();

public:
	ReorderingSSM(const label nnz,const label nRows, const label* row, const label* col);

	~ReorderingSSM();

	std::vector<label>* getVetexOrder();
	std::vector<label>* getEdgeOrder();
};

}

// Driver Code
#endif
