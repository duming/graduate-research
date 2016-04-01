#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string>
//#include <malloc.h>

#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <iterator>
#include <algorithm>

#define getmax(a,b) a>b?a:b
#define getmin(a,b) a>b?b:a




void PrintErrorAndQuit(char* sErrorString);


template <class A> void NewArray(A *** array, int Narray1, int Narray2)
{
  *array=new A* [Narray1];
  for(int i=0; i<Narray1; i++) *(*array+i)=new A [Narray2];
};

template <class A> void DeleteArray(A *** array, int Narray)
{
  for(int i=0; i<Narray; i++)
    if(*(*array+i)) delete [] *(*array+i);
  if(Narray) delete [] (*array);
  (*array)=NULL;
};

double dist(double x[3], double y[3]);

double dot(double *a, double *b);

void transform(double t[3], double u[3][3], double *x, double *x1);


void do_rotation(double **x, double **x1, int len, double t[3], double u[3][3]);





