/********************************************************************
	filename: 	PolyCubicSpline.h
	author:		hu zhijian
	created:	5:5:2010   16:52
	brief:	n元三次样条插值类的声明
*********************************************************************/

#ifndef NWPU_POLYCUBICSPLINE_H
#define NWPU_POLYCUBICSPLINE_H


#include <SL_dll.h>

#include <BaseInterp.h>
#include <CubicSpline.h>
#include <BiCubicSpline.h>
#include <vector>
#include <interp_struct.h>

namespace gslcpp
{
	class SL_DLL_API CPolyCubicSpline : public CBaseInterp 
	{
	private:
		CBaseInterp* method;                          // interpolation method
		double* csArray;                        // array for final cubic spline interpolation
		std::vector<CPolyCubicSpline*> pcs;               // array of PolyCubicSplineFasts for use with recursive step
		size_t dimOne;                                 // xArray dimension in a recursive step		
        double yValue;                           // returned interpolated value
		interp_table m_Table;
		size_t x_size; //要插值的数据x的个数
	public:

		CPolyCubicSpline(const interp_table& table);

		~CPolyCubicSpline();
	
		double RawInterp(const double* x, size_t jlo[]);

		double Interp(const double* x, const size_t& size);
	};

}


#endif