/********************************************************************
	filename: 	BaseInterp.cpp
	author:		huzhijian
	created:	5:5:2010   16:49
	brief:	implementation of the CBaseInterp class.
*********************************************************************/


#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <BaseInterp.h>
#include <gsl_errno.h>

namespace gslcpp
{
	//////////////////////////////////////////////////////////////////////
	// Construction/Destruction
	//////////////////////////////////////////////////////////////////////
	
	CBaseInterp::~CBaseInterp()
	{
		Free();
		
	}
	
	CBaseInterp::CBaseInterp(const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, 
			const double* y_Value, const size_t& y_size, const size_t& m, bool extro_flag) 
		:acc(NULL), jlo(NULL), xx_size(x_data_size), yy(y_Value), 
		 yy_size(y_size), dim(dimension), mm(m), extroInterp_flag(extro_flag)
	{
		
		Init(dimension, x_Value, x_data_size, y_size);
	}

	inline void CBaseInterp::Free()
	{
		
		if (jlo)
		{
			delete [] jlo;
			jlo = NULL;
		}
		if(acc)
		{
			delete [] acc;
			acc = NULL;
		}
		
	}
	size_t CBaseInterp::Locate(const double& x, const size_t& i)
	{
		size_t n = xx_size[i];
		const double* data = xx[i];

		size_t ju,jm,jl;
		
		if (n < 2 || mm < 2 || mm > n) 
		{
			gsl_error("size error!", __FILE__, __LINE__, GSL_EBADLEN);
		}
		
		bool ascnd=(data[n-1] >= data[0]); //True if ascending order of table, false otherwise.
		
		jl=0; //Initialize lower
		
		ju=n-1; //and upper limits.
		
		while (ju-jl > 1)  //二分法
		{ 
			//If we are not yet done,
			jm = (ju+jl) >> 1;/// compute a midpoint,
			
			if (x >= data[jm] == ascnd)
				jl=jm; //and replace either the lower limit
			else
				ju=jm; //or the upper limit, as appropriate.
			
		} //Repeat until the test condition is satisfied.
		
		acc[i].cor = (size_t)abs((int)jl- (int)(acc[i].cache)) > acc[i].dj ? 0 : 1; //Decide whether to use hunt or locate next time.
		//jl = __max( 0, __min( n-mm, jl-( (mm-2) >> 1) ) );
		acc[i].cache = jl;
		return __max( 0, __min( n-mm, jl-( (mm-2) >> 1) ) );
	}
	
	size_t CBaseInterp::Hunt(const double& x, const size_t& i)
	{
		size_t n = xx_size[i];
		const double* data = xx[i];

		size_t jl = acc[i].cache;
		size_t jm;
		size_t ju;
		size_t inc=1;

		if (n < 2 || mm < 2 || mm > n) 
		{
			gsl_error("size error!", __FILE__, __LINE__, GSL_EBADLEN);
		}

		bool ascnd=(data[n-1] >= data[0]);// True if ascending order of table, false otherwise.
		
		if (jl < 0 || jl > n-1)
		{ 
			//Input guess not useful. Go immediately to bisection.
			jl=0; 
			ju=n-1;
		} 
		else 
		{
			if ((x >= data[jl]) == ascnd) 
			{
				//Hunt up:
				for (;;) 
				{
					ju = jl + inc;
					if (ju >= n-1) { ju = n-1; break;} //Off end of table.
					else 
						if (x < data[ju] == ascnd) break;// Found bracket.
						else 
						{ 
							//Not done, so double the increment and try again.
							jl = ju;
							inc += inc;
						}
				}
			} 
			else 
			{
				// Hunt down:
				ju = jl;
				for (;;) 
				{
					if (jl <= inc) { jl = 0; break;} //Off end of table.
					jl = jl - inc;					
					if (x >= data[jl] == ascnd) break; //Found bracket.
					else 
					{ 
						//Not done, so double the increment and try again.
						ju = jl;
						inc += inc;	
					}
				}
			}
		}
		while (ju-jl > 1) 
		{ 
			//Hunt is done, so begin the final bisection phase:
			jm = (ju+jl) >> 1;
			if (x >= data[jm] == ascnd)
				jl=jm;
			else
				ju=jm;
		}
		acc[i].cor =  (size_t)abs((int)jl-(int)acc[i].cache) > acc[i].dj ? false : true;// Decide whether to use hunt or locate next time.
		acc[i].cache = jl; 
		return __max(0, __min( n-mm, jl- ( (mm-2) >> 1) ) );
	}

	double CBaseInterp::Interp(const double* x, const size_t& size)
	{
		//Given a value x, return an interpolated value, using data pointed to by xx and yy.
		if (dim != size)
		{
			gsl_error("interp函数参数x的维数不等于插值维数!", __FILE__, __LINE__, GSL_EBADLEN);
		}

		const double* data = NULL;
		for (size_t i=0; i<dim; i++)
		{
			data = xx[i];
			size_t n = xx_size[i];
			bool ascnd=(data[n-1] >= data[0]);

			if (!extroInterp_flag && (x[i]>data[n-1] == ascnd || x[i] < xx[i][0] == ascnd))
			{
				gsl_error("插值数据落在插值区间外边", __FILE__, __LINE__, GSL_EINVAL);
			}
			if (x[i] == xx[i][acc[i].cache]) //直接命中
			{
				jlo[i] = acc[i].cache;
			}
			else	//否则判断使用全区间二分法,还是局部二分法
			{
				jlo[i] = acc[i].cor ? Hunt(x[i], i) : Locate(x[i], i);
			}
		}
		
		
		return RawInterp(x, jlo);
		
	}
	
	void CBaseInterp::Init( const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, const size_t& y_size)
	{

		size_t total_num = 1;
		
		acc = new ACCEL[dimension]; //给加速器分配空间
		
		//存放自变量.
		for (size_t i=0; i<dimension; i++)
		{
			xx.push_back(x_Value[i]);
			
			//有可能溢出, fix me!
			total_num *= x_data_size[i];
			
			//加速器初始化
			acc[i].cache = 0;
			acc[i].cor = false;
			acc[i].dj = __min( 1, (int) pow( (double)x_data_size[i], 0.25) );
		}
		
		if (total_num != y_size )
		{
			gsl_error("the size of x and y do not conform", 
				__FILE__, __LINE__, GSL_EBADLEN );
		}
		
		jlo = new size_t[dim];
	}
}
