/********************************************************************
	filename: 	CubicSpline.cpp
	author:		huzhijian
	created:	5:5:2010   16:51
	brief:	一元 三次样条插值类
*********************************************************************/

#include <stdlib.h>
#include <CubicSpline.h>
#include <MathConst.h>
#include <gsl_errno.h>

namespace gslcpp
{
	
	CCubicSpline::CCubicSpline(const interp_table& table) 
		: CBaseInterp(1, table.x_Array, table.x_Array_Num, table.y_Array, table.y_Array_Num, 2, table.extro_interp_flag), //一元两点
		derivCalculated(false), nPoints(*(table.x_Array_Num)), x(table.x_Array[0]), y(table.y_Array)
	{
		if(nPoints<3)
		{
			gsl_error("A minimum of three data points is needed ", __FILE__, __LINE__, GSL_EINVAL);
		}
		d2ydx2.resize(nPoints);

	    CalcDeriv();
	}

	CCubicSpline::~CCubicSpline()
	{
		Free();
	}

	void CCubicSpline::Free()
	{
	
	}

	void CCubicSpline::CalcDeriv()
	{
		double	p=0.0,qn=0.0,sig=0.0,un=0.0;
		
		double* u = new double[nPoints];

		//const double* x = xx[0]->GetDataPtr();
		//const double* y = yy.GetDataPtr();

		d2ydx2[0] = u[0] = 0.0;

		

		for(size_t i=1;i<=this->nPoints-2;i++)
		{
			sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
			p=sig*d2ydx2[i-1]+2.0;
			d2ydx2[i]=(sig-1.0)/p;
			u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
			u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		}

		

		qn=un=0.0;
		
		this->d2ydx2[nPoints-1]=(un-qn*u[nPoints-2])/(qn*d2ydx2[nPoints-2]+1.0);
		
		for(size_t k=nPoints-2, j=0; j<=nPoints-2;k--, j++)
		{
			d2ydx2[k]=d2ydx2[k]*d2ydx2[k+1]+u[k];
		}
		delete []  u;
	    this->derivCalculated = true;
	}

	double CCubicSpline::RawInterp(const double* xArray, size_t jlo[])
	{
		double x_to_intep = xArray[0];	//给定插值点
		size_t index = jlo[0];	//插值点所在的位置
		double y_res;
		double h=0.0,b=0.0,a=0.0;
		size_t khi=index+1;
	    size_t klo=index;
		h = x[khi]-x[klo];
		
		if (h == 0.0)
		{
			printf("Two values of x are identical: point x[%u]=x[%u] = %g", index, index+1, x[index]);
			abort();
		}
		else
		{
			a=(x[khi] - x_to_intep)/h;
			b=(x_to_intep - x[klo])/h;
			y_res = a*y[klo]+b*y[khi]+((a*a*a-a)*d2ydx2[klo]+(b*b*b-b)*d2ydx2[khi])*(h*h)/6.0;
			
		}
	    return y_res;
	}
}