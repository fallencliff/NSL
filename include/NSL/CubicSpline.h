/********************************************************************
	filename: 	CubicSpline.h
	author:		hu zhijian
	created:	5:5:2010   16:51
	brief:	һԪ ����������ֵ��
*********************************************************************/

#ifndef NSL_CUBICSPLINE_H__
#define NSL_CUBICSPLINE_H__


#include <NSL.H>
#include <BaseInterp.h>
#include <interp_struct.h>
#include <vector>


namespace gslcpp
{
	class NSL_EXPORT CCubicSpline : public CBaseInterp 
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



#endif // NSL_CUBICSPLINE_H__