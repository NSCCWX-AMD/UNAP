#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//- test for sizeof

using namespace std;

void main()
{
	double* p = new double[10];

	int size1 = sizeof(*p);
	int size2 = sizeof(double);
	int size3 = sizeof(p);

	cout << "size of *p is " << size1 << endl;
	cout << "size of double is " << size2 << endl;
	cout << "size of pointer is " << size3 << endl;


	class A
	{
	private:
		int i_;
		double v_;

	public:
		A(int i, double v)
		:
			i_(i),
			v_(v)
		{}
	};

	A* pA = new A(2, -0.2);
	A aa(2, -0.2);
	A& aa2 = *pA;

	size1 = sizeof(pA);
	size2 = sizeof(aa);
	size3 = sizeof(*pA);
	int size4 = sizeof(aa2);

	cout << "size of pointer is " << size1 << endl;
	cout << "size of object is " << size2 << endl;
	cout << "size of *pointer is " << size3 << endl;
	cout << "size of reference is " << size4 << endl;

}