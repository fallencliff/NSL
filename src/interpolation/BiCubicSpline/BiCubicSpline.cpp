/********************************************************************
	filename: 	BiCubicSpline.cpp
	author:		huzhijian
	created:	5:5:2010   16:50
	brief:	二元 三次样条插值类
*********************************************************************/

#include <stdio.h> 
#include <BiCubicSpline.h>
#include <CubicSpline.h>
#include <MathConst.h>
#include <gsl_errno.h>

namespace gslcpp
{

	CBiCubicSpline::~CBiCubicSpline()
	{
 		if (yRow_interp)
 		{
 			delete [] yRow_interp;
			yRow_interp= NULL;
 		}
		
		for(int i=0; i<nPoints; i++)
		{
			if (srp[i])
			{
				delete srp[i];
				srp[i] = NULL;
			}
		}
		
	}

	CBiCubicSpline::CBiCubicSpline(const interp_table& table)
		: CBaseInterp(2, table.x_Array, table.x_Array_Num, table.y_Array, table.y_Array_Num, 2, table.extro_interp_flag),
			nPoints(table.x_Array_Num[0]), x(table.x_Array[0]), yRow_interp(NULL)

	{
 		yRow_interp = new double[nPoints];
		srp.resize(nPoints);
		global_table = table;

		interp_table row_intrp = table;
		row_intrp.dimension =1;
		row_intrp.x_Array = &(table.x_Array[1]);
		row_intrp.x_Array_Num = &table.x_Array_Num[1];

		for(size_t i=0; i<nPoints; i++)
		{
			row_intrp.y_Array = &table.y_Array[i*table.x_Array_Num[1]];

			row_intrp.y_Array_Num = table.x_Array_Num[1];

			srp[i] = new CCubicSpline(row_intrp); //construct CCubicSpline on each row.
	    }
	}

	double CBiCubicSpline::RawInterp(const double* xArray, size_t jlo[])
	{
		for (size_t i=0; i<nPoints; i++)
		{
			yRow_interp[i] = srp[i]->Interp(&xArray[1], 1);//Interpolate on each row.
		}

		interp_table col_interp = global_table;
		col_interp.dimension =1;
		col_interp.y_Array = yRow_interp;
		col_interp.y_Array_Num = nPoints;
		CCubicSpline scol(col_interp); //Construct the column spline, and evaluate it
		
		return scol.Interp(&xArray[0], 1);
	}
	
	double CBiCubicSpline::Interp(const double* x, const size_t& size)
	{
		return RawInterp(x, 0); //重写此函数是为了减少下标寻找运算，0为占位符
	}
}