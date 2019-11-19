#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//- test for deleting class pointer
//- and for creating an object in construct function

using namespace std;

class A
{
private:
	int    a1_;
	double a2_;

public:
	A()
	:
		a1_(100),
		a2_(1e-06)
	{}

	A
	(
		int a1,
		double a2
	)
	:
		a1_(a1),
		a2_(a2)
	{}

	~A()
	{
		cout << "I am A" << endl;
	}

	int a1() const
	{
		return a1_;
	}

	double a2() const
	{
		return a2_;
	}
};

class B
{
	bool deleteAPtr_;
	A *APtr_;

public:
	B()
	:
		deleteAPtr_(false)
	{
		A *a1 = new A;
		APtr_ = a1;
		deleteAPtr_ = true;
	}

	B(A &a)
	:
		deleteAPtr_(false),
		APtr_(&a)
	{}

	~B()
	{
		if(deleteAPtr_)
		{
			delete APtr_;
			APtr_ = NULL;
		}
		cout << "I am B" << endl;
	}

	void dispA()
	{
		cout << "In a, a1 = " << APtr_->a1() << ", a2 = " << APtr_->a2() << endl;
	}
};

int main()
{
	A aa(5, 0.88);

	B bb(aa);

	bb.dispA();

	bb.~B();

	cout << "After destruction, a1 = " << aa.a1() << ", a2 = " << aa.a2() << endl;

	cout << "I am finished" << endl;

	return 0;
}

