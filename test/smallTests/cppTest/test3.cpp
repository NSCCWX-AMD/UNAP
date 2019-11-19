#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//- test for virtual function but return different object

using namespace std;

class A
{
private:

public:

	virtual void display() const
	{
		cout << "I am A here." << endl;
		cout << "so the answer is NO!" << endl;
	}

};

class B
:
	public A
{
private:

public:

	virtual void display() const
	{
		cout << "I am B here." << endl;
		cout << "so the answer is YES!" << endl;
	}

};

class C
{
private:

public:

	virtual A& show() = 0;

};


class D
:
	public C
{
private:
	B* bb_;

public:

	D()
	:
		bb_(NULL)
	{}

	virtual B& show()
	{
		bb_ = new B;
		return *bb_;
	}

	virtual ~D()
	{
		if(bb_)
		{
			delete bb_;
			bb_ = NULL;
		}
	}
};


int main()
{
	D dd;

	A &aa = dd.show();

	aa.display();

	return 0;
}



