/********************************************************************
	filename: 	PolyCubicSpline.cpp
	author:		huzhijian
	created:	5:5:2010   16:52
	brief:	n元三次样条插值类
*********************************************************************/

#include <CubicSpline.h>
#include <BiCubicSpline.h>
#include <PolyCubicSpline.h>

namespace gslcpp
{

	CPolyCubicSpline::CPolyCubicSpline(const interp_table& table )//: CBaseInterp(x, y, dimension, 2)
		: CBaseInterp(table.dimension, table.x_Array, table.x_Array_Num, table.y_Array, table.y_Array_Num, 2, table.extro_interp_flag),
		  m_Table(table), x_size(0), csArray(NULL), method(NULL)
			//nPoints(table.x_Array_Num[0]), x(table.x_Array[0]), yRow_interp(NULL)
	{
		double yValue = 0.0;

		size_t dim = GetNumOfDimension();
		switch(dim)
		{
		case 1:    
			{
				// 若为一维插值
				
				this->method = new CCubicSpline(m_Table);							
				break;
			}

		case 2: 
			{
				// 若为二维插值
				this->method = new CBiCubicSpline(m_Table);							
				break;
			}	

		default:   
			// 大于二维时，递归调用CPolyCubicSpline
			dimOne = table.x_Array_Num[0];
			csArray = new double [dimOne];

			size_t block_num=1;

			for (size_t n=1; n<dim; n++)
			{
				block_num *= table.x_Array_Num[n];
			}
			const double* data = table.y_Array;
			const size_t col_num = table.x_Array_Num[table.dimension-1];
			interp_table temp = m_Table;
			temp.dimension = m_Table.dimension - 1;
			for(size_t i=0; i<dimOne; i++)
			{		
				temp.x_Array = &(table.x_Array[1]);
				temp.x_Array_Num = &(table.x_Array_Num[1]);
				temp.y_Array = table.y_Array+ i* block_num;
				temp.y_Array_Num = block_num;
				
				
				CPolyCubicSpline* p = new CPolyCubicSpline(temp);
				pcs.push_back(p);
				
			}
         }
	}


	CPolyCubicSpline::~CPolyCubicSpline()
	{
		//清理已申请的内存
		if (csArray)
		{
			delete csArray;
			csArray = NULL;
		}
		
		size_t size = pcs.size();
		for(size_t i=0; i<size; i++)
		{
			if (pcs[i])
			{
				delete pcs[i];
				pcs[i] = NULL;
			}
		}
		if (method)
		{
			delete method;
			method = NULL;
		}
		
	}

	double CPolyCubicSpline::RawInterp(const double* x, size_t jlo[])
	{	
		switch(GetNumOfDimension())
		{
		//case 0:     throw new IllegalArgumentException("data array must have at least one dimension");
		case 1:    
			// If it is one dimensional perform simple cubic spline
			this->yValue = ((CCubicSpline*)(this->method))->Interp(x, 1);

			break;


		case 2:     
			// If it is two dimensional perform bicubic spline
			this->yValue = ((CBiCubicSpline*)(this->method))->Interp(x, 2);

			break;


		default:   
			// If it is greater than two dimensional, recursively call PolyCubicSplineFast
			//CVectorView newCoord = xArray.Subvector(1, dim-1);

			for(int i=0; i<this->dimOne; i++)
			{
				csArray[i] = pcs[i]->Interp(&x[1], x_size-1);
			}
			
			// Perform simple cubic spline on the array of above returned interpolates

			interp_table temp = m_Table;

			temp.dimension = 1;
			temp.y_Array = csArray;
			temp.y_Array_Num = dimOne;

			CCubicSpline ncs(temp);

			this->yValue = ncs.Interp(x, 1);
			
		}
		
        return this->yValue;

	}
	
	double CPolyCubicSpline::Interp(const double* x, const size_t& size)
	{
		x_size = size; 
		return RawInterp(x, 0);
	}

}

