#include <stdio.h>
#include <stdlib.h>
#include <iostream>

//- test for PtrList

using namespace std;

template<typename T>
class PtrList
{
private:

	//- number of T
	int size_;

	//- T* pointer
	T** PtrList_;

public:

	PtrList();

	PtrList(int size);

	~PtrList()
	{
		for(int i=0; i<size_; i++)
		{
			cout << " i = " << i << ", size = " << size_ << endl;
			if(PtrList_[i])
			{
				delete PtrList_[i];
				PtrList_[i] = NULL;
			}
		}

		delete [] PtrList_;
		PtrList_ = NULL;
		size_ = 0;
	}


	int size() const
	{
		return size_;
	}

	//- return element const reference
    inline const T* operator[](const int i) const
    {
    	return PtrList_[i];
    }

    //- return element pointer reference
    T& operator[](const int i)
    {
    	return *PtrList_[i];
    }

    void removeLevel(const int i);

    void setLevel(const int i, T& obj);

    void SET_size(const int i);
};


class A
{
private:

	int    a1_;
	double a2_;
	int* jPtr_;

public:
	A()
	:
		a1_(100),
		a2_(1e-06),
		jPtr_(NULL)
	{
		jPtr_ = new int[a1_];
		for(int i=0; i<a1_; i++)
		{
			jPtr_[i] = i;
		}
	}

	A
	(
		int a1,
		double a2
	)
	:
		a1_(a1),
		a2_(a2),
		jPtr_(NULL)
	{
		jPtr_ = new int[a1_];
		for(int i=0; i<a1_; i++)
		{
			jPtr_[i] = i;
		}
	}

	~A()
	{
		cout << "I am A.a1 = " << a1_ << endl;
		delete []  jPtr_;
	}

	int &a1()
	{
		return a1_;
	}

	double &a2()
	{
		return a2_;
	}
};


template<typename T>
PtrList<T>::PtrList()
:
	size_(50),
	PtrList_(NULL)
{
	PtrList_ = new T* [size_];
	for(int i=0; i<size_; i++)
	{
		PtrList_[i] = NULL;
	}
}


template<typename T>
PtrList<T>::PtrList(int size)
:
	size_(size),
	PtrList_(NULL)
{
	PtrList_ = new T* [size_];
	for(int i=0; i<size_; i++)
	{
		PtrList_[i] = NULL;
	}
}

template<typename T>
void PtrList<T>::removeLevel(const int i)
{
	if(i < size_)
	{
		// (*PtrList_[i]).~T();
		delete PtrList_[i];
		for(int j=i; j<size_-1; j++)
		{
			PtrList_[j] = PtrList_[j+1];
		}
		PtrList_[size_ - 1] = NULL;
        size_--;
	}
	else
	{
		cout << "removeLevel error" << endl;
	}
}

template<typename T>
void PtrList<T>::setLevel(const int i, T& obj)
{
	if(i < size_)
	{
		PtrList_[i] = &obj;
	}
	else
	{
		cout << "setLevel error" << endl;
	}
}

template<typename T>
void PtrList<T>::SET_size(const int newSize)
{
	int oldSize = size();
	if(newSize == 0)
	{
		this->~PtrList();
	}
	else if(oldSize == newSize)
	{
		// nothing to do
	}
	else
	{
		T** PtrListOld = PtrList_;
		PtrList_ = new T* [newSize];

		for(int i=0; i<newSize; i++)
		{

			if(i < oldSize)
			{
				PtrList_[i] = PtrListOld[i];
			}
			else
			{
				PtrList_[i] = NULL;
			}
		}
		size_ = newSize;
		delete [] PtrListOld;
	}
}


int main()
{
	A* a1 = new A(1, 0.1);
	A* a2 = new A(2, 0.2);
	A* a3 = new A(3, 0.3);
	A* a4 = new A(4, 0.4);

	PtrList<A> aa;
	cout << "aa size is " << aa.size() << endl;

	// PtrList<A> aa(4);
	aa.setLevel(0, *a1);
	aa.setLevel(1, *a2);
	aa.setLevel(2, *a3);
	aa.setLevel(3, *a4);

	aa.SET_size(0);
	cout << "aa size is " << aa.size() << endl;

	cout << "Before: " << endl;
	for(int i=0; i<aa.size(); i++)
	{
		cout << "At i = " << i << ", a1 is " << aa[i].a1() << ", a2 is " << aa[i].a2() << endl;
	}

	aa.removeLevel(1);

	cout << "After: " << endl;
	for(int i=0; i<aa.size(); i++)
	{
		cout << "At i = " << i << ", a1 is " << aa[i].a1() << ", a2 is " << aa[i].a2() << endl;
	}

	return 0;
}
