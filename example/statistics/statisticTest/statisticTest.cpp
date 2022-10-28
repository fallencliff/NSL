/********************************************************************
	filename: 	statisticTest.cpp
	author:		hu zhijian
	created:	14:5:2010   9:25
	brief:	Here is a basic example of how to use the statistical class
*********************************************************************/


#include <stdio.h>
#include <Statis.h>

using gslcpp::CStatis;
int main(int argc, char* argv[])
{
	double data[5] = {17.2, 18.1, 16.5, 18.3, 12.6};
	double mean, variance, largest, smallest;

	CStatis st(data, 5, 1);

	mean = st.get_mean(); //¾ùÖµ
	variance = st.variance(); //·½²î
	largest = st.max();
	smallest = st.min();

	printf ("The dataset is %g, %g, %g, %g, %g\n", data[0], data[1], data[2], data[3], data[4]);

	printf ("The sample mean is %g\n", mean);
	printf ("The estimated variance is %g\n", variance);
	printf ("The largest value is %g\n", largest);
	printf ("The smallest value is %g\n", smallest);

	return 0;
}

