/********************************************************************
	filename: 	PolyLagrangeInterp.h
	author:		hu zhijian
	created:	5:5:2010   16:53
	brief:	多元拉格朗日插值类
*********************************************************************/
#ifndef NWPU_POLYLAGRANGEINTERP_H__
#define NWPU_POLYLAGRANGEINTERP_H__

#include <SL_dll.h>

#include <BaseInterp.h>
#include <interp_struct.h>
#include <vector>

namespace gslcpp
{
	class SL_DLL_API CPolyLagrangeInterp : public CBaseInterp  
	{
	public:

		CPolyLagrangeInterp(const interp_table& interp_data);

		virtual ~CPolyLagrangeInterp();

	private:	
		double RawInterp(const double* x, size_t jlo[]);

		double Mul(const double* x, double t, size_t m, size_t i); //连乘

		double Add(size_t next); //连加

		std::vector<size_t> global_array;	//每个连加当前的下标

		const double* x_to_interp; //给定点

		std::vector<size_t> index_base; //用来确定函数值的位置

	};
}


#endif // POLYLAGRANGEINTERP_H__
