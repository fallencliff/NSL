// vectorTest.cpp : Defines the entry point for the console application.
//
/********************************************************************
	filename: 	Load.cpp
	author:		hu zhijian
	created:	13:5:2010   15:42
	brief:	an example show how to read data from file to vector
*********************************************************************/

#include <iostream>
#include <valarray>
#include <MathConst.h>
#include <Vector.h>

using std::valarray;
using namespace gslcpp;
using std::cout;
using std::endl;

template<typename T>
void showVector(const Vector<T>& v)
{
	for (size_t i=0; i<v.size(); i++)
	{
		cout << v.cref(i) << " ";
	}
	cout << endl;
}
int main(int argc, char* argv[])
{
	
	
	Vector<double> v(3);

	v<<1,2,3;

	v = v*3.0;

	showVector(v);

//	valarray<double> val = v.toVal();
	cout << absMax(v);

	v.zero();

	showVector(v);

	Vector<double> v1;

	cout<< CMath::max(1,2);

	return 0;
}

