// vectorTest.cpp : Defines the entry point for the console application.
//
/********************************************************************
	filename: 	vectorTest.cpp
	author:		hu zhijian
	created:	13:5:2010   9:09
	brief:	This program shows how to allocate, initialize and read from a vector
*********************************************************************/

#include <stdio.h>
#include <Vector.h>

using namespace gslcpp;
int main(int argc, char* argv[])
{
	
	
	CVector v(3);

	int i;
	for (i = 0; i < 3; i++)
	{
		v[i] = 1.23 + i;
	}

	for (i = 0; i < 3; i++)
	{
		printf ("v_%d = %g\n", i,v[i]);
	}

	return 0;
}

