#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <typeinfo>

#include "mpi.h"
#include "unap.hpp"

using namespace std;

template <typename T>
void reduce(T &v)
{
  T vLocal = v;
  CommData mytype;
  if (typeid(T) == typeid(int))
  {
    cout << "COMM_INT" << endl;
    mytype = COMM_INT;
  }
  else if (typeid(T) == typeid(double))
  {
    cout << "COMM_DOUBLE" << endl;
    mytype = COMM_DOUBLE;
  }
  else if (typeid(T) == typeid(long int))
  {
    cout << "COMM_LONG" << endl;
    mytype = COMM_LONG;
  }
  else if (typeid(T) == typeid(float))
  {
    cout << "COMM_FLOAT" << endl;
    mytype = COMM_FLOAT;
  }
  UNAP::unapMPI::unapCommunicator().allreduce(
      "allreduce", &vLocal, &v, 1, mytype, COMM_SUM);
  UNAP::unapMPI::unapCommunicator().finishTask("allreduce");
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

  if (myid == 0)
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

  if (myid == 0)
  {
    cout << "mpi: " << endl;
    cout << "COMM_INT size = " << sizeof(COMM_INT) << endl;
    cout << "COMM_LONG size = " << sizeof(COMM_LONG) << endl;
    cout << "COMM_FLOAT size = " << sizeof(COMM_FLOAT) << endl;
    cout << "COMM_DOUBLE size = " << sizeof(COMM_DOUBLE) << endl;

    cout << "compiler: " << endl;
    cout << "INT size = " << sizeof(int) << endl;
    cout << "LONG size = " << sizeof(long int) << endl;
    cout << "FLOAT size = " << sizeof(float) << endl;
    cout << "DOUBLE size = " << sizeof(double) << endl;

    cout << "a = " << a << ", b = " << b << ", c = " << c << ", d = " << d
         << endl;
  }

  MPI_Finalize();

  return 0;
}
