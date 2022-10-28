/********************************************************************
	filename: 	MathConst.cpp
	author:		hu zhijian
	created:	30:5:2010   17:15
	brief:	
*********************************************************************/


#include <MathConst.h>

namespace gslcpp
{
	namespace CMath
	{	
		double 	fdiv (const double x, const double y)
		{
			return x / y;
		}
		
		inline double nan (void)
		{
			return fdiv (0.0, 0.0);
		}
		
		inline double posinf (void)
		{
			return fdiv (+1.0, 0.0);
		}
		
		inline double neginf (void)
		{
			return fdiv (-1.0, 0.0);
		}
		
		
		inline int is_nan (double x)
		{
			return _isnan(x);
		}
		
		inline int is_inf (double x)
		{
			int fpc = _fpclass(x);
			
			if (fpc == _FPCLASS_PINF)
				return +1;
			else if (fpc == _FPCLASS_NINF)
				return -1;
			else 
				return 0;
		}
		
		inline int is_finite ( double x )
		{
			return _finite(x);
		}
		
		inline double log1p ( double x )
		{
			volatile double y;
			y = 1 + x;
			return log(y) - ((y-1)-x)/y ;  /* cancels errors with IEEE arithmetic */
		}
		
		inline double expm1 ( double x )
		{
			/* FIXME: this should be improved */
			
			if (fabs(x) < LN2)
			{
				/* Compute the taylor series S = x + (1/2!) x^2 + (1/3!) x^3 + ... */
				
				double i = 1.0;
				double sum = x;
				double term = x / 1.0;
				
				do
				{
					i++ ;
					term *= x/i;
					sum += term;
				}
				while (fabs(term) > fabs(sum) * DBL_EPSILON) ;
				
				return sum ;
			}
			else
			{
				return exp(x) - 1;
			}
		}
		
		inline double hypot ( double x, double y )
		{
			double xabs = fabs(x) ;
			double yabs = fabs(y) ;
			double min, max;
			
			if (xabs < yabs) {
				min = xabs ;
				max = yabs ;
			} else {
				min = yabs ;
				max = xabs ;
			}
			
			if (min == 0) 
			{
				return max ;
			}
			
			{
				double u = min / max ;
				return max * sqrt (1 + u * u) ;
			}
		}
		
		inline double acosh ( double x )
		{
			if (x > 1.0 / DBL_EPSILON)
			{
				return log (x) + LN2;
			}
			else if (x > 2)
			{
				return log (2 * x - 1 / (sqrt (x * x - 1) + x));
			}
			else if (x > 1)
			{
				double t = x - 1;
				return log1p (t + sqrt (2 * t + t * t));
			}
			else if (x == 1)
			{
				return 0;
			}
			else
			{
				return NAN;
			}
		}
		
		inline double asinh ( double x )
		{
			double a = fabs (x);
			double s = (x < 0) ? -1 : 1;
			
			if (a > 1 / SQRT_DBL_EPSILON)
			{
				return s * (log (a) + LN2);
			}
			else if (a > 2)
			{
				return s * log (2 * a + 1 / (a + sqrt (a * a + 1)));
			}
			else if (a > SQRT_DBL_EPSILON)
			{
				double a2 = a * a;
				return s * log1p (a + a2 / (1 + sqrt (1 + a2)));
			}
			else
			{
				return x;
			}
		}
		
		inline double atanh ( double x )
		{
			double a = fabs (x);
			double s = (x < 0) ? -1 : 1;
			
			if (a > 1)
			{
				return NAN;
			}
			else if (a == 1)
			{
				return (x < 0) ? NEGINF : POSINF;
			}
			else if (a >= 0.5)
			{
				return s * 0.5 * log1p (2 * a / (1 - a));
			}
			else if (a > DBL_EPSILON)
			{
				return s * 0.5 * log1p (2 * a + 2 * a * a / (1 - a));
			}
			else
			{
				return x;
			}
		}
		
		inline double ldexp ( double x, int e )
		{
			double p2 = pow(2.0, (double)e);
			return x * p2;
		}
		inline double frexp ( double x, int* e )
		{
			if (x == 0.0)
			{
				*e = 0;
				return 0.0;
			}
			else
			{
				double ex = ceil (log (fabs (x)) / LN2);
				int ei = (int) ex;
				double f = ldexp (x, -ei);
				
				while (fabs (f) >= 1.0)
				{
					ei++;
					f /= 2.0;
				}
				
				while (fabs (f) < 0.5)
				{
					ei--;
					f *= 2.0;
				}
				
				*e = ei;
				return f;
			}
		}
		
		inline double pow_int ( double x, int n )
		{
			double value = 1.0;
			
			if(n < 0) {
				x = 1.0/x;
				n = -n;
			}
			
			/* repeated squaring method 
			* returns 0.0^0 = 1.0, so continuous in x
			*/
			do {
				if(n & 1) value *= x;  /* for n odd */
				n >>= 1;
				x *= x;
			} while (n);
			
			return value;
		}
		
		
		inline double pow_2 ( double x )
		{
			return x*x;
		}
		
		inline double pow_3 ( double x )
		{
			return x*x*x;
		}
		
		inline double pow_4 ( double x )
		{
			double x2 = x*x;  
			return x2*x2; 
		}
		
		inline double pow_5 ( double x )
		{
			double x2 = x*x;   
			return x2*x2*x;
		}
		
		inline double pow_6 ( double x )
		{
			double x2 = x*x;   
			return x2*x2*x2;
		}
		
		inline double pow_7 ( double x )
		{
			double x3 = x*x*x; 
			return x3*x3*x;
		}
		
		inline double pow_8 ( double x )
		{
			double x2 = x*x;   
			double x4 = x2*x2; 
			return x4*x4; 
		}
		
		inline double pow_9 ( double x )
		{
			double x3 = x*x*x; 
			return x3*x3*x3;
		}
		inline int sign ( float x )
		{
			return SIGN(x);
		}
		
		inline int sign ( double x )
		{
			return SIGN(x);
		}
		
		inline int sign ( int x )
		{
			return SIGN(x);
		}
		
		inline int sign ( long x )
		{
			return SIGN(x);
		}
		
		inline bool is_odd ( int x )
		{
			return IS_ODD(x);
		}
		
		inline bool is_even ( int x )
		{
			return IS_EVEN(x);
		}
		
		inline int max ( int x, int y )
		{
			return MAX (x, y);
		}
		
		inline double max ( double x, double y )
		{
			return MAX (x, y);
		}
		
		inline long double max ( long double x, long double y )
		{
			return MAX (x, y);
		}
		
		inline int min ( int x, int y )
		{
			return MIN (x, y);
		}
		
		inline double min ( double x, double y )
		{
			return MIN (x, y);
		}
		
		inline long double min ( long double x, long double y )
		{
			return MIN (x, y);
		}
		
		inline int fcmp ( double x1, double x2, double epsilon )
		{
			
			int exponent;
			double delta, difference;
			
			/* Find exponent of largest absolute value */
			
			{
				double max = (fabs (x1) > fabs (x2)) ? x1 : x2;
				
				frexp (max, &exponent);
			}
			
			/* Form a neighborhood of size  2 * delta */
			
			delta = ldexp (epsilon, exponent);
			
			difference = x1 - x2;
			
			if (difference > delta)       /* x1 > x2 */
			{
				return 1;
			}
			else if (difference < -delta) /* x1 < x2 */
			{
				return -1;
			}
			else                          /* -delta <= difference <= delta */
			{
				return 0;                 /* x1 ~=~ x2 */
			}
		}
	}
}
