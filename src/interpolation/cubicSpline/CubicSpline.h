/********************************************************************
	filename: 	CubicSpline.h
	author:		hu zhijian
	created:	5:5:2010   16:51
	brief:	һԪ ����������ֵ��
*********************************************************************/

#ifndef	NWPU_CUBICSPLINE_H
#define NWPU_CUBICSPLINE_H

#include <SL_dll.h>
#include <BaseInterp.h>
#include <interp_struct.h>
#include <vector>


namespace gslcpp
{
	class SL_DLL_API CCubicSpline : public CBaseInterp 
	{

		friend class CBiCubicSpline;
	private:
		size_t nPoints;                 // x�ĸ���	
		std::vector<double> d2ydx2;					//���׵���
		bool derivCalculated;          // = true �����Ѿ������	
		const double* x;						// x in f(x)
		const double* y;						// y=f(x) tabulated function
	private:
		void Free();

	protected:

	public:
		CCubicSpline(const interp_table& table);

		virtual ~CCubicSpline();
	
		//  Calculates the second derivatives
		void CalcDeriv();

		double RawInterp(const double* xArray, size_t jlo[]);
		
	};

}

#endif