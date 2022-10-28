// vectorTest.cpp : Defines the entry point for the console application.
//
/********************************************************************
	filename: 	save.cpp
	author:		hu zhijian
	created:	13:5:2010   15:40
	brief:	The program shows how to write a vector to a file.
*********************************************************************/

#include <stdio.h>
#include <vector.h>

using namespace gslcpp;

int main(int argc, char* argv[])
{
	
	CVector v(100);

	size_t i;

	for (i = 0; i < 100; i++)
	{
		v[i] = 1.23 + i;
	}

	v.Save("test.dat", "%g");

	return 0;
}

