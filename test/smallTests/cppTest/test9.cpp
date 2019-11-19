#include <memory>   //- shared_ptr, unique_ptr, weak_ptr
#include <iostream>

using namespace std;

class A
{
private:
	int i_;
	int* iPtr_;


public:
	A(int i)
	:
		i_(i),
		iPtr_(NULL)
	{
		iPtr_ = new int[i_];
	}

	~A()
	{
		cout << "Here is the destructor of A." << endl;
		if(iPtr_)
		{
			cout << "here ok" << endl;
			delete []iPtr_;
			iPtr_ = NULL;
		}
	}
};


class B
{
private:
	shared_ptr<A> Aptr_;

public:
	B(shared_ptr<A> APtr)
	:
		Aptr_(APtr)
	{

	}
	~B()
	{
		cout << "Here is the destructor of B." << endl;
	}
};



int main(int argc, char const *argv[])
{
	shared_ptr<A> APtr(new A(2));
	B bb(APtr);
	// delete APtr;
	APtr->~A();
	cout << "middle" << endl;

	return 0;
}






