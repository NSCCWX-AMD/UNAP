// C++ program for Implementation of
// Reverse Cuthill Mckee Algorithm
#ifndef cooRCM_HPP_
#define cooRCM_HPP_
#include <bits/stdc++.h>

namespace COORCM
{
class ReorderingSSM {
private:
	const int _nnz;    // number of non-zeros
	const int _nRows;  // number of rows
	const int* _row;   // array of rows
	const int* _col;   // array of columns

	static std::vector<int> _globalDegree;
	std::vector<int>* _vetexOrderPtr;
	std::vector<int>* _edgeOrderPtr;

	// int* _rowsOffset;

	void degreeGenerator();
    static bool compareDegree(int i, int j);

	void CuthillMckee();

	// Implementation of reverse Cuthill-Mckee algorithm
	void ReverseCuthillMckee();

public:
	ReorderingSSM(const int nnz,const int nRows, const int* row, const int* col);

	~ReorderingSSM();

	std::vector<int>* getVetexOrder();
	std::vector<int>* getEdgeOrder();
};

}

// Driver Code
#endif
