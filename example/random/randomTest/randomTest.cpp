/********************************************************************
	filename: 	randomTest.cpp
	author:		hu zhijian
	created:	13:5:2010   21:22
	brief:	The following program demonstrates the use of a random number generator to produce
			uniform random numbers in the range [0.0, 1.0), and [0, n]
*********************************************************************/


#include <stdio.h>
#include <RandGenerator.h>
#include <MathConst.h>
using namespace gslcpp;

// int main(int argc, char* argv[])
// {
// 
// 	size_t i=0;
// 
// 	CRandom def; //use default generator
// 
// 	FILE* fp = fopen("test.dat", "w");
// 
// 	if (!fp)
// 	{
// 		return -1;
// 	}
// 
// 	for (i = 0; i < 10000; i++)
// 	{
// 		double u = def.RandUniform(); //use default seed
// 
// 		//fprintf (fp, "%f\t%f\n", u, u);
// 		fprintf (fp, "%f\n", u);
// 	}
// 
// 	fclose(fp);
// 
// 	return 0;
// }

int main(int argc, char* argv[])
{
	
	size_t i=0;
	
	CRandom rng(RNG_RANLUX); 
	
	FILE* fp = fopen("test.dat", "w");
	
	if (!fp)
	{
		return -1;
	}
	
	printf("%s\n", rng.Name());
	for (i = 0; i < 10000; i++)
	{
		double u = rng.RandUniform(); //use default seed
		
		//fprintf (fp, "%f\t%f\n", u, u);
		fprintf (fp, "%f\n", u);
	}
	
	fclose(fp);

	return 0;
}