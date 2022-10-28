/********************************************************************
	filename: 	BiCubicSpline.h
	author:		hu zhijian
	created:	5:5:2010   16:49
	brief:	��Ԫ ����������ֵ��
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
		const double* x;					//��ŵ�һά����ָ��
		double* yRow_interp;						//ÿһ�в�ֵ�Ľ��
		std::vector<CCubicSpline*> srp;		//���һά��ֵ����ָ��
		interp_table global_table;
	public:
		virtual ~CBiCubicSpline();

		CBiCubicSpline(const interp_table& table);
		double RawInterp(const double* xArray, size_t jlo[]);
		double Interp(const double* x, const size_t& size);
	};

}

#endif