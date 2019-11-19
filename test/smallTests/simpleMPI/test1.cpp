#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <typeinfo>

using namespace std;

template<typename T>
void reduce(T& v)
{
	T vLocal = v;
	MPI_Datatype mytype;
	if(typeid(T) == typeid(int))
	{
		cout << "MPI_INT" << endl;
		mytype = MPI_INT;
	}
	else if(typeid(T) == typeid(double))
	{
		cout << "MPI_DOUBLE" << endl;
		mytype = MPI_DOUBLE;
	}
	else if(typeid(T) == typeid(long int))
	{
		cout << "MPI_LONG" << endl;
		mytype = MPI_LONG;
	}
	else if(typeid(T) == typeid(float))
	{
		cout << "MPI_FLOAT" << endl;
		mytype = MPI_FLOAT;
	}
	MPI_Allreduce(&vLocal, &v, 1, mytype, MPI_SUM, MPI_COMM_WORLD);
}


int main(int argc, char *argv[])
{
	int myid, num_procs;
   	MPI_Init(&argc, &argv);
   	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
   	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

   	int a;
   	double b;
   	long int c;
   	float d;

   	if(myid == 0)
   	{
   		a = 1;
   		b = 0.1;
   		c = 100;
   		d = 1.1;
   	}
   	else
   	{
   		a = 2;
   		b = 0.2;
   		c = 200;
   		d = 2.2;
   	}

   	// if(num_procs > 1)
   	{
   		reduce(a);
   		reduce(b);
   		reduce(c);
   		reduce(d);
   	}

   	if(myid == 0)
   	{
   		cout << "mpi: " << endl;
   		cout << "MPI_INT size = " << sizeof(MPI_INT) << endl;
   		cout << "MPI_LONG size = " << sizeof(MPI_LONG) << endl;
   		cout << "MPI_FLOAT size = " << sizeof(MPI_FLOAT) << endl;
   		cout << "MPI_DOUBLE size = " << sizeof(MPI_DOUBLE) << endl;

   		cout << "compiler: " << endl;
   		cout << "INT size = " << sizeof(int) << endl;
   		cout << "LONG size = " << sizeof(long int) << endl;
   		cout << "FLOAT size = " << sizeof(float) << endl;
   		cout << "DOUBLE size = " << sizeof(double) << endl;

   		cout << "a = " << a << ", b = " << b
   			 << ", c = " << c << ", d = " << d << endl;
   	}

   	MPI_Finalize();


	return 0;
}



