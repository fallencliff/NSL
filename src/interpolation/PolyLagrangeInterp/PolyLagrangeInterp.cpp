/********************************************************************
	filename: 	PolyLagrangeInterp.cpp
	author:		huzhijian
	created:	5:5:2010   16:55
	brief:	implementation of the CPolyLagrangeInterp class.
*********************************************************************/


#include <PolyLagrangeInterp.h>
#include <math.h>
#include <string.h>

namespace gslcpp
{
	//////////////////////////////////////////////////////////////////////
	// Construction/Destruction
	//////////////////////////////////////////////////////////////////////
	
	CPolyLagrangeInterp::~CPolyLagrangeInterp()
	{

	}
	
	CPolyLagrangeInterp::CPolyLagrangeInterp(const interp_table& interp_data)
	: CBaseInterp(interp_data.dimension, interp_data.x_Array, interp_data.x_Array_Num, interp_data.y_Array,
				  interp_data.y_Array_Num, interp_data.Lagrange_point_Num, interp_data.extro_interp_flag)
	{
		size_t dim = interp_data.dimension;

		index_base.resize(dim);

		if (dim >=3)
		{
			index_base[dim-2] = 1;
			for (size_t i=dim-3,j=0; j<dim-2; i--, j++)
			{
				index_base[i] = index_base[i+1] * interp_data.x_Array_Num[i+1];
			}
		}
		else
		{
			index_base[0] = dim==2 ? 1: 0;
		}
		
		
	}

	double CPolyLagrangeInterp::RawInterp(const double* x, size_t jlo[])
	{

		double res = 0.0;

		x_to_interp = x;

		const size_t dim = GetNumOfDimension();

		for (size_t index=0; index<dim; index++)
		{
			global_array.push_back(jlo[index]);//[index] = jlo[index];
		}

		return Add(dim-1); //dim个连加
	}

	//局部连乘函数
	double CPolyLagrangeInterp::Mul(const double* x, double t, size_t m, size_t i)
	{
		double res = 1.0;
		const size_t mm = Get_mm();
		for(size_t j=m; j<=m+mm-1; j++)
		{
			
			if(j == i)
			{
				
				continue;
				
			}
			res = res * (t - x[j]) / (x[i] - x[j]);
			
		}
		
		return res;
	}

	double CPolyLagrangeInterp::Add(size_t add_next)
	{
		const size_t* jlo = Get_jlo_Addr();
		const size_t dim = GetNumOfDimension();
		const size_t mm = Get_mm();
		const double* y_addr = Get_Y_Addr();

		if (0 == add_next)
		{
			global_array[add_next] = jlo[add_next];
			double add_res = 0.0;
			do 
			{	
				double mul_res = 1.0;

				for (size_t i=0; i<dim; i++) //dim个连乘
				{
					mul_res *= Mul(Get_X_Addr(i), x_to_interp[i], jlo[i], global_array[i]);
				}

				size_t col_index = global_array[dim-1];
				size_t row_index = 0;
				
				for (size_t j=0; j<dim-1; j++)
				{					
					row_index += global_array[j]* index_base[j];
				}
				
				mul_res *= y_addr[row_index*Get_X_Size(dim-1) + col_index];//此处要乘函数值

				global_array[add_next]++;

				add_res += mul_res;

			} while (global_array[add_next] <= jlo[add_next]+mm-1);
			
			return add_res; 
			
		}

		global_array[add_next] = jlo[add_next]; //初始化连加的初始下标

		double sum = 0.0;

		do 
		{
			sum += Add(add_next-1);

			global_array[add_next]++;

		} while (global_array[add_next] <= jlo[add_next]+mm-1); //每个连加做mm次
		//jlo[add_next]+mm-1为连加的结束下标

		return sum;
	}
}
