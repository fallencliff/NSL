// mathTest.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <RandGenerator.h>
#include <MathConst.h>

//using gslcpp::CMath;
using namespace gslcpp;

int main(int argc, char* argv[])
{
	printf("pi = %.15f\n", CMath::PI);
	printf("e = %.15f\n", CMath::E);

	CRandom rng;
	printf("%g\n", rng.RandUniform());
	printf("%g\n", rng.RandUniform());
	printf("%g\n", rng.Rand());
	return 0;
}

