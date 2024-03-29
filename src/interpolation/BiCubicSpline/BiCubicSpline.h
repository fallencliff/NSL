/********************************************************************
	filename: 	BiCubicSpline.h
	author:		hu zhijian
	created:	5:5:2010   16:49
	brief:	二元 三次样条插值类
*********************************************************************/


#ifndef	NWPU_BICUBICSPLINE_H
#define NWPU_BICUBICSPLINE_H

#include <SL_dll.h>

#include <BaseInterp.h>
#include <CubicSpline.h>
#include <interp_struct.h>
#include <vector>



namespace gslcpp
{
	class SL_DLL_API CBiCubicSpline : public CBaseInterp 
	{
	private:	
		size_t nPoints;   	                            // no. of x1 tabulated points
		const double* x;					//存放第一维数据指针
		double* yRow_interp;						//每一行插值的结果
		std::vector<CCubicSpline*> srp;		//存放一维插值对象指针
		interp_table global_table;
	public:
		virtual ~CBiCubicSpline();

		CBiCubicSpline(const interp_table& table);
		double RawInterp(const double* xArray, size_t jlo[]);
		double Interp(const double* x, const size_t& size);
	};

}

#endif