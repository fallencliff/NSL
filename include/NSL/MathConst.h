/********************************************************************
filename: 	MathConst.h
author:		hu zhijian
created:	12:5:2010   21:51
brief:	
*********************************************************************/

#ifndef NSL_MATHCONST_H__
#define NSL_MATHCONST_H__
#include <float.h>
#include <NSL.h>
#include <math.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define IS_ODD(n)  ((n) & 1)
#define IS_EVEN(n) (!(IS_ODD(n)))
#define SIGN(x)    ((x) >= 0.0 ? 1 : -1)


namespace gslcpp
{
	#ifndef NSL_CMath_H__
	#define NSL_CMath_H__

	namespace  CMath
	{
#define   E 2.71828182845904523536028747135
		
#define LOG2E  1.44269504088896340735992468100      /* log_2 (e) */
		
#define LOG10E  0.43429448190325182765112891892      /* log_10 (e) */
		
#define SQRT2  1.41421356237309504880168872421     /* sqrt(2) */
		
#define SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
		
#define SQRT3  1.73205080756887729352744634151      /* sqrt(3) */
		
#define PI  3.14159265358979323846264338328      /* pi */
		
#define PI_2  1.57079632679489661923132169164      /* pi/2 */
		
#define PI_4  0.78539816339744830966156608458      /* pi/4 */
		
#define SQRTPI  1.77245385090551602729816748334     /* sqrt(pi) */
		
#define D2_SQRTPI  1.12837916709551257389615890312      /* 2/sqrt(pi) */
		
#define D1_PI  0.31830988618379067153776752675     /* 1/pi */
		
#define D2_PI  0.63661977236758134307553505349      /* 2/pi */
		
#define LN10  2.30258509299404568401799145468     /* ln(10) */
		
#define LN2  0.69314718055994530941723212146      /* ln(2) */
		
#define LNPI  1.14472988584940017414342735135     /* ln(pi) */
		
#define EULER  0.57721566490153286060651209008      /* Euler static constant */
#define DOUBLE_EPSILON  2.2204460492503131e-16
#define SQRT_DBL_EPSILON   1.4901161193847656e-08
#define ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define LOG_DBL_EPSILON   -3.6043653389117154e+01
		
#define DOUBLE_MIN        2.2250738585072014e-308
#define SQRT_DBL_MIN   1.4916681462400413e-154
#define ROOT3_DBL_MIN  2.8126442852362996e-103
#define ROOT4_DBL_MIN  1.2213386697554620e-77
#define ROOT5_DBL_MIN  2.9476022969691763e-62
#define ROOT6_DBL_MIN  5.3034368905798218e-52
#define LOG_DBL_MIN   -7.0839641853226408e+02
		
#define DOUBLE_MAX        1.7976931348623157e+308
#define SQRT_DBL_MAX   1.3407807929942596e+154
#define ROOT3_DBL_MAX  5.6438030941222897e+102
#define ROOT4_DBL_MAX  1.1579208923731620e+77
#define ROOT5_DBL_MAX  4.4765466227572707e+61
#define ROOT6_DBL_MAX  2.3756689782295612e+51
#define LOG_DBL_MAX    7.0978271289338397e+02
		
		double posinf (void);
		double neginf (void);
		double nan (void);
#define POSINF  posinf()
		
#define NEGINF  neginf()
		
#define NAN  nan()
		
			/** 
			 * @brief 如果x不是数字返回1.	
			 *
			 * @org: gsl_isnan
			 */
		NSL_EXPORT	 int is_nan (double x);

			/** 
			 * @brief x是正无穷返回+1, 负无穷返回-1, 否则返回0.	
			 *
			 * @org: gsl_isinf
			 */
		NSL_EXPORT	 int is_inf (double x);

			/** 
			 * @brief 如果x是实数返回1, 无穷大或者不是数字则返回0.	
			 *
			 * @org: gsl_finite
			 */
		NSL_EXPORT	 int is_finite (double x);

			/** 
			 * @brief $log(1+x)$.	
			 *
			 * @org: gsl_log1p
			 */
		NSL_EXPORT	 double log1p (double x);

			/** 
			 * @brief $exp(x)-1$.	
			 *
			 * @org: gsl_expm1
			 */
		NSL_EXPORT	 double expm1 (double x);

			/** 
			 * @brief $sqrt{x^2 + y^2}$.	
			 *
			 * @org: gsl_hypot
			 */
		NSL_EXPORT	 double hypot (double x, double y);

			/** 
			 * @brief $arccosh(x)$.	
			 *
			 * @org: gsl_acosh
			 */
		NSL_EXPORT	 double acosh (double x);

			/** 
			 * @brief $arcsinh(x)$.	
			 *
			 * @org: gsl_asinh
			 */
		NSL_EXPORT	 double asinh (double x);

			/** 
			 * @brief $arctanh(x)$.	
			 *
			 * @org: gsl_atanh
			 */
		NSL_EXPORT	 double atanh (double x);

			/** 
			 * @brief $x * 2^e$.	
			 *
			 * @org: gsl_ldexp
			 */
		NSL_EXPORT	 double ldexp (double x, int e);

			/** 
			 * @brief 此函数把x分解为一个系数f和一个指数e, 例如 $x = f * 2^e$ 并且 $0.5 <= f < 1$.	
			 *
			 * @org: gsl_frexp
			 */
		NSL_EXPORT	 double frexp (double x, int* e);

			/** 
			 * @brief $x^n$.	
			 *
			 * @org: gsl_pow_int
			 */
		NSL_EXPORT	 double pow_int (double x, int n);

			/** 
			 * @brief $x^2$.	
			 *
			 * @org: gsl_pow_2
			 */
		NSL_EXPORT	 double pow_2 (double x);

			/** 
			 * @brief $x^3$.	
			 *
			 * @org: gsl_pow_3
			 */
		NSL_EXPORT	 double pow_3 (double x);

			/** 
			 * @brief $x^4$.	
			 *
			 * @org: gsl_pow_4
			 */
		NSL_EXPORT	 double pow_4 (double x);

			/** 
			 * @brief $x^5$.	
			 *
			 * @org: gsl_pow_5
			 */
		NSL_EXPORT	 double pow_5 (double x);

			/** 
			 * @brief $x^6$.	
			 *
			 * @org: gsl_pow_6
			 */
		NSL_EXPORT	 double pow_6 (double x);

			/** 
			 * @brief $x^7$.	
			 *
			 * @org: gsl_pow_7
			 */
		NSL_EXPORT	 double pow_7 (double x);

			/** 
			 * @brief $x^8$.	
			 *
			 * @org: gsl_pow_8
			 */
		NSL_EXPORT	 double pow_8 (double x);

			/** 
			 * @brief $x^9$.	
			 *
			 * @org: gsl_pow_9
			 */
		NSL_EXPORT	 double pow_9 (double x);

			/** 
			 * @brief 返回x的符号. 定义如下((x) >= 0 ? 1 : -1).  注意在这种定义下零的符号是正的.	
			 *
			 * @org: SIGN
			 */
		NSL_EXPORT	 int sign (float x);

			/** 
			 * @brief 返回x的符号. 定义如下((x) >= 0 ? 1 : -1).  注意在这种定义下零的符号是正的.	
			 *
			 * @org: SIGN
			 */
		NSL_EXPORT	 int sign (double x);

			/** 
			 * @brief 返回x的符号. 定义如下((x) >= 0 ? 1 : -1).  注意在这种定义下零的符号是正的.	
			 *
			 * @org: SIGN
			 */
		NSL_EXPORT	 int sign (int x);

			/** 
			 * @brief 返回x的符号. 定义如下((x) >= 0 ? 1 : -1).  注意在这种定义下零的符号是正的.	
			 *
			 * @org: SIGN
			 */
		NSL_EXPORT	 int sign (long x);

			/** 
			 * @brief 若x为奇数返回true.	
			 *
			 * @org: IS_ODD
			 */
		NSL_EXPORT	 bool is_odd (int x);

			/** 
			 * @brief 若x为偶数返回true.	
			 *
			 * @org: IS_EVEN
			 */
		NSL_EXPORT	 bool is_even (int x);

			/** 
			 * @brief 返回x和y的最大值.		
			 *
			 * @org: GSL_MAX_DBL
			 */
		NSL_EXPORT	 double max (double x, double y);

			/** 
			 * @brief 返回x和y的最大值.		
			 *
			 * @org: GSL_MAX_LDBL
			 */
		NSL_EXPORT	 long double max (long double x, long double y);

			/** 
			 * @brief 返回x和y的最大值.		
			 *
			 * @org: GSL_MAX_INT
			 */
		NSL_EXPORT	 int max (int x, int y);

			/** 
			 * @brief 返回x和y的最大值.	
			 *
			 * @org: GSL_MIN_DBL
			 */
		NSL_EXPORT	 double min (double x, double y);

			/** 
			 * @brief 返回x和y的最小值.		
			 *
			 * @org: GSL_MIN_LDBL
			 */
		NSL_EXPORT	 long double min (long double x, long double y);

			/** 
			 * @brief 返回x和y的最小值.		
			 *
			 * @org: GSL_MIN_INT
			 */
		NSL_EXPORT	 int min (int x, int y);

			/** 
			 * @brief 双精度浮点数x和y的比较.		
			 *
			 * @return 如果x和y在epsilon允许的精度下大约相等则返回0.
			 * @		如果 x < y, 返回 -1.
			 * @		如果 x > y, 返回 +1. 
			 *
			 * @org: gsl_fcmp
			 */
		NSL_EXPORT	 int fcmp (double x1, double x2, double epsilon);
	}

	#endif // NSL_CMath_H__
}




#endif // NSL_MATHCONST_H__