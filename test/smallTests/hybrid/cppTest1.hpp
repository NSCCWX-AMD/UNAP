#ifndef CPPTEST1_HPP
#define CPPTEST1_HPP

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

class A
{
protected:
	int a1_;
	double a2_;

public:
	A(int i, double v);
	~A()
	{}

	virtual void display() const = 0;
};

class B
:
	public A
{
private:

public:
	B(int i, double v);
	~B()
	{}

	virtual void display() const;
};

#endif
