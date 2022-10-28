/********************************************************************
	filename: 	MatrixTest.cpp
	author:		hu zhijian
	created:	13:5:2010   15:48
	brief:	The program below shows how to allocate, initialize and read from a matrix
*********************************************************************/


#include <stdio.h>
#include <Matrix.h>

using gslcpp::CMatrix;
using gslcpp::CMatrixView;

void PrintMatrix(const CMatrix & m)
{
	for (size_t i=0; i<m.Size1(); i++)
	{
		for (size_t j=0; j<m.Size2(); j++)
		{
			printf("%g ", m[i][j]);
		}
		
		printf("\n");
	}

	printf("\n");
}

int main()
{

	CMatrix mv(10, 10);

	mv.Load("a.txt");

	//PrintMatrix(mv);

	CMatrix res = mv.Inv();

	PrintMatrix(res);

	
	
}

/*
int main(int argc, char* argv[])
{

	CMatrix matrix(10, 3);

	size_t i, j;
	for (i=0; i<matrix.Size1(); i++)
	{
		for (j=0; j<matrix.Size2(); j++)
		{
			matrix[i][j] = 0.23 + 100*i + j;
		}

	}
	for (i=0; i<matrix.Size1(); i++)
	{
		for (j=0; j<matrix.Size2(); j++)
		{
			printf("m(%d,%d) = %g\n", i, j, matrix[i][j]);
		}
		
	}
	
	double data[6] = {1,2,3,4,5,6};
	
	CMatrix m1(data, 1, 6); //1*6 矩阵
	CMatrix m2(data, 6, 1); //6*1
	return 0 ;
}
*/
/************************************************************************/
/* The next program shows how to write and read a matrix to or from a file 
/************************************************************************/

/*
int main (int argc, char* argv[])
{
	size_t i, j, k1 = 0, k2 =0;

	size_t size1 = 50;
	size_t size2 = 50;
	CMatrix ma(size1, size2);

	CMatrix mb(size1, size2), mc(size1, size2);

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			ma[i][j] = 0.23 + i + j;
		}
		
	}

	ma.Save("testText.dat"); //文本方式输出
	ma.SaveBinary("testBinary.dat"); //二进制方式输出
	
	mc.Load("testText.dat"); //文本方式输入
	mb.LoadBinary("testBinary.dat"); //二进制方式输入
	

	for (i = 0; i < size1; i++)
	{
		for (j = 0; j < size2; j++)
		{
			if (ma[i][j] != mb[i][j])
				k1++;

			if (ma[i][j] != mc[i][j])
				k2++;
		}
	}
	
	printf ("differences = %d (should be zero)\n", k1);
	printf ("differences = %d (might not be zero)\n", k2);

	return (k1 > 0);
}

*/