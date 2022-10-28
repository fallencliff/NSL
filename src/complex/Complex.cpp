/********************************************************************
	filename: 	Complex.cpp
	author:		huzhijian
	created:	5:5:2010   16:51
	brief:	复数类
*********************************************************************/

#include <math.h>
#include <MathConst.h>
#include <Complex.h>

namespace gslcpp
{
#define GSL_REAL(z)     ((z).dat[0])
#define GSL_IMAG(z)     ((z).dat[1])
#define GSL_COMPLEX_P(zp) ((zp)->dat)
#define GSL_COMPLEX_P_REAL(zp)  ((zp)->dat[0])
#define GSL_COMPLEX_P_IMAG(zp)  ((zp)->dat[1])
#define GSL_COMPLEX_EQ(z1,z2) (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))
	
#define GSL_SET_COMPLEX(zp,x,y) do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#define GSL_SET_REAL(zp,x) do {(zp)->dat[0]=(x);} while(0)
#define GSL_SET_IMAG(zp,y) do {(zp)->dat[1]=(y);} while(0)

	CComplex::CComplex()
	{		
		GSL_SET_COMPLEX (&gsldata, 0.0, 0.0);	
	}

	CComplex::CComplex(const double r, const double i)
	{
		GSL_SET_COMPLEX (&gsldata, r, i);
	}

	CComplex::CComplex(const CComplex& c)
	{		
		gsldata = c.gsldata;	
	}

	CComplex::CComplex(gsl_complex z) : gsldata(z) 
	{ 
	}

	CComplex CComplex::rect(double x, double y)
	{
		return CComplex(x, y);
	}

	CComplex CComplex::polar(double r, double theta)
	{
		return CComplex(r * ::cos (theta), r * ::sin (theta));
	}

	double CComplex::real() const
	{
		return GSL_REAL(gsldata);
	}

	double CComplex::imag() const
	{
		return GSL_IMAG(gsldata);
	}

	void CComplex::set(double x, double y)
	{
		GSL_SET_COMPLEX(&gsldata, x, y);
	}

	void CComplex::set_real(double a)
	{
		GSL_SET_REAL(&gsldata, a);
	}

	void CComplex::set_imag(double a)
	{
		GSL_SET_IMAG(&gsldata, a);
	}

	// 基础函数
	double CComplex::arg() const
	{
		double x = GSL_REAL (gsldata);
		double y = GSL_IMAG (gsldata);
		
		if (x == 0.0 && y == 0.0)
		{
			return 0;
		}
		
		return ::atan2 (y, x);
	}

	double CComplex::abs() const
	{
		return ::hypot (GSL_REAL (gsldata), GSL_IMAG (gsldata));
	}

	double CComplex::abs2() const
	{
		double x = GSL_REAL (gsldata);
		double y = GSL_IMAG (gsldata);
		
		return (x * x + y * y);
	}

	double CComplex::logabs() const
	{
		double xabs = fabs (GSL_REAL (gsldata));
		double yabs = fabs (GSL_IMAG (gsldata));
		double max, u;
		
		if (xabs >= yabs)
		{
			max = xabs;
			u = yabs / xabs;
		}
		else
		{
			max = yabs;
			u = xabs / yabs;
		}
		
		/* Handle underflow when u is close to 0 */
		
		return ::log (max) + 0.5 * CMath::log1p (u * u);
	}
	
	CComplex CComplex::pow(CComplex& c) const
	{
		gsl_complex z;
		
		if (GSL_REAL (gsldata) == 0 && GSL_IMAG (gsldata) == 0.0)
		{
			GSL_SET_COMPLEX (&z, 0.0, 0.0);
		}
		else
		{
			double logr = logabs();
			double theta = arg ();
			
			double br = GSL_REAL (c.gsldata), bi = GSL_IMAG (c.gsldata);
			
			double rho = ::exp (logr * br - bi * theta);
			double beta = theta * br + bi * logr;
			
			GSL_SET_COMPLEX (&z, rho * ::cos (beta), rho * ::sin (beta));
		}
	
		return CComplex(z);
	}
	
	CComplex CComplex::pow(double x) const
	{
		gsl_complex z;
		
		if (GSL_REAL (gsldata) == 0 && GSL_IMAG (gsldata) == 0)
		{
			GSL_SET_COMPLEX (&z, 0, 0);
		}
		else
		{
			double logr = logabs ();
			double theta = arg ();
			double rho = ::exp (logr * x);
			double beta = theta * x;
			GSL_SET_COMPLEX (&z, rho * ::cos (beta), rho * ::sin (beta));
		}
		return CComplex(z);
	}

	CComplex CComplex::sqrt()
	{
		gsl_complex z;
		
		if (GSL_REAL (gsldata) == 0.0 && GSL_IMAG (gsldata) == 0.0)
		{
			GSL_SET_COMPLEX (&z, 0, 0);
		}
		else
		{
			double x = fabs (GSL_REAL (gsldata));
			double y = fabs (GSL_IMAG (gsldata));
			double w;
			
			if (x >= y)
			{
				double t = y / x;
				w = ::sqrt (x) * ::sqrt (0.5 * (1.0 + ::sqrt (1.0 + t * t)));
			}
			else
			{
				double t = x / y;
				w = ::sqrt (y) * ::sqrt (0.5 * (t + ::sqrt (1.0 + t * t)));
			}
			
			if (GSL_REAL (gsldata) >= 0.0)
			{
				double ai = GSL_IMAG (gsldata);
				GSL_SET_COMPLEX (&z, w, ai / (2.0 * w));
			}
			else
			{
				double ai = GSL_IMAG (gsldata);
				double vi = (ai >= 0) ? w : -w;
				GSL_SET_COMPLEX (&z, ai / (2.0 * vi), vi);
			}
		}	 
		return CComplex(z);
	}

	CComplex CComplex::sqrt(double x)
	{
		gsl_complex z;
		
		if (x >= 0)
		{
			GSL_SET_COMPLEX (&z, ::sqrt (x), 0.0);
		}
		else
		{
			GSL_SET_COMPLEX (&z, 0.0, ::sqrt (-x));
		}
		return CComplex(z);
	}

	CComplex CComplex::log10() const
	{
		return log()*(1/::log(10.0));
	}

	CComplex CComplex::log() const
	{
		double logr = logabs ();
		double theta = arg ();
		
		return CComplex(logr, theta);
	}

	CComplex CComplex::log_b(const CComplex& b) const
	{
		return log()/b.log();
	}

	CComplex CComplex::exp() const
	{
		double rho = ::exp (GSL_REAL (gsldata));
		double theta = GSL_IMAG (gsldata);
		
		return CComplex(rho * ::cos (theta), rho * ::sin (theta));
	}
	// 重载操作符
	CComplex CComplex::operator+(const CComplex& c) const
	{	
		return CComplex(real()+c.real(), imag()+c.imag());
	}

	CComplex CComplex::operator-(const CComplex& c) const
	{
		return CComplex(real()-c.real(), imag()-c.imag());
	}

	CComplex CComplex::operator*(const CComplex& c) const
	{
		double ar = GSL_REAL (gsldata), ai = GSL_IMAG (gsldata);
		double br = GSL_REAL (c.gsldata), bi = GSL_IMAG (c.gsldata);
		
		return CComplex(ar * br - ai * bi, ar * bi + ai * br);
	}

	CComplex CComplex::operator/(const CComplex& c) const
	{
		double ar = GSL_REAL (gsldata), ai = GSL_IMAG (gsldata);
		double br = GSL_REAL (c.gsldata), bi = GSL_IMAG (c.gsldata);
		
		double s = 1.0 / c.abs();
		
		double sbr = s * br;
		double sbi = s * bi;
		
		double zr = (ar * sbr + ai * sbi) * s;
		double zi = (ai * sbr - ar * sbi) * s;
  
		return CComplex(zr, zi);
	}

	CComplex CComplex::operator+(double x) const
	{
		return CComplex(real()+x, imag());
	}

	CComplex CComplex::operator-(double x) const
	{
		return CComplex(real()-x, imag());
	}

	CComplex CComplex::operator*(double x) const
	{
		return CComplex(real()*x, imag()*x);
	}

	CComplex CComplex::operator/(double x) const
	{
		return CComplex(real()/x, imag()/x);
	}

	CComplex CComplex::operator=(double x)
	{
		this->set(x, 0.0);
		
		return *this;
	}

	CComplex CComplex::operator=(const CComplex& c)
	{
		this->set(c.real(), c.imag());

		return *this;
	}

	int CComplex::operator==(const CComplex& c) const
	{
		int res = 1;
		
		res = res && (CMath::fcmp(this->real(), c.real(),10.0*DOUBLE_EPSILON) == 0) && 
				 (CMath::fcmp(this->imag(), c.imag(),10.0*DOUBLE_EPSILON) == 0);
		   
		return res;
		
	}

	int CComplex::operator==(double x) const
	{
		int res = 1;
		res = res && (CMath::fcmp(this->real(), x,10.0*DOUBLE_EPSILON) == 0)&& 
			(CMath::fcmp(this->imag(), 0.0,10.0*DOUBLE_EPSILON) == 0);
		
		return res;
		
	}

	CComplex CComplex::add_imag(double y) const
	{
		return CComplex(real(), imag()+y);
	}

	CComplex CComplex::sub_imag(double y) const
	{
		return CComplex(real(), imag()-y);
	}

	CComplex CComplex::mul_imag(double y) const
	{
		return CComplex(-y*imag(), y*real());
	}

	CComplex CComplex::div_imag(double y) const
	{
		return CComplex(imag()/y, real()/(-y));
	}

	CComplex CComplex::conjugate() const
	{
		return CComplex(real(), -1*imag());
	}

	CComplex CComplex::inverse() const
	{
		double s = 1.0 / abs();
		return CComplex(real()*s*s, -1*imag()*s*s);
	}

	CComplex CComplex::operator-(void) const
	{
		return CComplex(-1*real(), -1*imag());
	}

	// 三角函数
	CComplex CComplex::sin() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;
		
		if (I == 0.0) 
		{
			/* avoid returing negative zero (-0.0) for the imaginary part  */
			
			GSL_SET_COMPLEX (&z, ::sin (R), 0.0);  
		} 
		else 
		{
			GSL_SET_COMPLEX (&z, ::sin (R) * ::cosh (I), ::cos (R) * ::sinh (I));
		}
		
		return CComplex(z);
	}

	CComplex CComplex::cos() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;
		
		if (I == 0.0) 
		{
			/* avoid returing negative zero (-0.0) for the imaginary part  */
			
			GSL_SET_COMPLEX (&z, ::cos (R), 0.0);  
		} 
		else 
		{
			GSL_SET_COMPLEX (&z, ::cos (R) * ::cosh (I), ::sin (R) * ::sinh (-I));
		}
		return CComplex(z);
	}

	CComplex CComplex::tan() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;
		
		if (fabs (I) < 1)
		{
			double D = ::pow (::cos (R), 2.0) + ::pow (::sinh (I), 2.0);
			
			GSL_SET_COMPLEX (&z, 0.5 * ::sin (2 * R) / D, 0.5 * ::sinh (2 * I) / D);
		}
		else
		{
			double u = ::exp (-I);
			double C = 2 * u / (1 - ::pow (u, 2.0));
			double D = 1 + ::pow (::cos (R), 2.0) * ::pow (C, 2.0);
			
			double S = ::pow (C, 2.0);
			double T = 1.0 / ::tanh (I);
			
			GSL_SET_COMPLEX (&z, 0.5 * ::sin (2 * R) * S / D, T / D);
		}

		return CComplex(z);
	}

	CComplex CComplex::sec() const
	{
		return cos().inverse();
	}

	CComplex CComplex::csc() const
	{
		return sin().inverse();
	}

	CComplex CComplex::cot() const
	{
		return tan().inverse();
	}

	// 反三角函数
	CComplex CComplex::arcsin() const
	{

		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		gsl_complex z;
		
		if (I == 0)
		{
			z = arcsin(R).gsldata;
		}
		else
		{
			double x = fabs (R), y = fabs (I);
			double r = hypot (x + 1, y), s = hypot (x - 1, y);
			double A = 0.5 * (r + s);
			double B = x / A;
			double y2 = y * y;
			
			double real, imag;
			
			const double A_crossover = 1.5, B_crossover = 0.6417;
			
			if (B <= B_crossover)
			{
				real = asin (B);
			}
			else
			{
				if (x <= 1)
				{
					double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
					real = ::atan (x / ::sqrt (D));
				}
				else
				{
					double Apx = A + x;
					double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
					real = ::atan (x / (y * ::sqrt (D)));
				}
			}
			
			if (A <= A_crossover)
			{
				double Am1;
				
				if (x < 1)
				{
					Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
				}
				else
				{
					Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
				}
				
				imag = CMath::log1p (Am1 + ::sqrt (Am1 * (A + 1)));
			}
			else
			{
				imag = ::log (A + ::sqrt (A * A - 1));
			}
			
			GSL_SET_COMPLEX (&z, (R >= 0) ? real : -real, (I >= 0) ? imag : -imag);
		}
				
		return CComplex(z);
	}

	CComplex CComplex::arcsin(double x)
	{
		gsl_complex z;
		
		if (fabs (x) <= 1.0)
		{
			GSL_SET_COMPLEX (&z, asin (x), 0.0);
		}
		else
		{
			if (x < 0.0)
			{
				GSL_SET_COMPLEX (&z, -1*PI_2, CMath::acosh (-x));
			}
			else
			{
				GSL_SET_COMPLEX (&z, PI_2, -1*CMath::acosh (x));
			}
		}

		return CComplex(z);
	}

	CComplex CComplex::arccos() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		gsl_complex z;
		
		if (I == 0)
		{
			z = arccos(R).gsldata;
		}
		else
		{
			double x = fabs (R), y = fabs (I);
			double r = hypot (x + 1, y), s = hypot (x - 1, y);
			double A = 0.5 * (r + s);
			double B = x / A;
			double y2 = y * y;
			
			double real, imag;
			
			const double A_crossover = 1.5, B_crossover = 0.6417;
			
			if (B <= B_crossover)
			{
				real = acos (B);
			}
			else
			{
				if (x <= 1)
				{
					double D = 0.5 * (A + x) * (y2 / (r + x + 1) + (s + (1 - x)));
					real = ::atan (::sqrt (D) / x);
				}
				else
				{
					double Apx = A + x;
					double D = 0.5 * (Apx / (r + x + 1) + Apx / (s + (x - 1)));
					real = ::atan ((y * ::sqrt (D)) / x);
				}
			}
			
			if (A <= A_crossover)
			{
				double Am1;
				
				if (x < 1)
				{
					Am1 = 0.5 * (y2 / (r + (x + 1)) + y2 / (s + (1 - x)));
				}
				else
				{
					Am1 = 0.5 * (y2 / (r + (x + 1)) + (s + (x - 1)));
				}
				
				imag = CMath::log1p (Am1 + ::sqrt (Am1 * (A + 1)));
			}
			else
			{
				imag = ::log (A + ::sqrt (A * A - 1));
			}
			
			GSL_SET_COMPLEX (&z, (R >= 0) ? real : PI - real, (I >= 0) ? -imag : imag);
		}
		
		return CComplex(z);
	}

	CComplex CComplex::arccos(double x)
	{
		gsl_complex z;
		
		if (fabs (x) <= 1.0)
		{
			GSL_SET_COMPLEX (&z, ::acos (x), 0);
		}
		else
		{
			if (x < 0.0)
			{
				GSL_SET_COMPLEX (&z, PI, -1*CMath::acosh (-x));
			}
			else
			{
				GSL_SET_COMPLEX (&z, 0,  CMath::acosh (x));
			}
		}
		return CComplex(z);
	}
	CComplex CComplex::arctan() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		gsl_complex z;
		
		if (I == 0)
		{
			GSL_SET_COMPLEX (&z, atan (R), 0);
		}
		else
		{
		/* FIXME: This is a naive implementation which does not fully
		take into account cancellation errors, overflow, underflow
			etc.  It would benefit from the Hull et al treatment. */
			
			double r = hypot (R, I);
			
			double imag;
			
			double u = 2 * I / (1 + r * r);
			
			/* FIXME: the following cross-over should be optimized but 0.1
			seems to work ok */
			
			if (fabs (u) < 0.1)
			{
				imag = 0.25 * ( CMath::log1p (u) -  CMath::log1p (-u));
			}
			else
			{
				double A = ::hypot (R, I + 1);
				double B = ::hypot (R, I - 1);
				imag = 0.5 * ::log (A / B);
			}
			
			if (R == 0)
			{
				if (I > 1)
				{
					GSL_SET_COMPLEX (&z,  PI_2, imag);
				}
				else if (I < -1)
				{
					GSL_SET_COMPLEX (&z, -1* PI_2, imag);
				}
				else
				{
					GSL_SET_COMPLEX (&z, 0, imag);
				};
			}
			else
			{
				GSL_SET_COMPLEX (&z, 0.5 * atan2 (2 * R, ((1 + r) * (1 - r))), imag);
			}
		}

		return CComplex(z);
	}

	CComplex CComplex::arcsec() const
	{
		return  inverse().arccos();
	}

	CComplex CComplex::arcsec(double x)
	{
		gsl_complex z;
		
		if (x <= -1.0 || x >= 1.0)
		{
			GSL_SET_COMPLEX (&z, acos (1 / x), 0.0);
		}
		else
		{
			if (x >= 0.0)
			{
				GSL_SET_COMPLEX (&z, 0,  CMath::acosh (1 / x));
			}
			else
			{
				GSL_SET_COMPLEX (&z, PI, -1* CMath::acosh (-1 / x));
			}
		}

		return CComplex(z);
	}
	CComplex CComplex::arccsc() const
	{
		return inverse().arcsin();
	}

	CComplex CComplex::arccsc(double x)
	{
		gsl_complex z;
		
		if (x <= -1.0 || x >= 1.0)
		{
			GSL_SET_COMPLEX (&z, asin (1 / x), 0.0);
		}
		else
		{
			if (x >= 0.0)
			{
				GSL_SET_COMPLEX (&z,  PI_2, -1* CMath::acosh (1 / x));
			}
			else
			{
				GSL_SET_COMPLEX (&z, -1* PI_2,  CMath::acosh (-1 / x));
			}
		}
		return CComplex(z);
	}

	CComplex CComplex::arccot() const
	{
		gsl_complex z;
		
		if (GSL_REAL (gsldata) == 0.0 && GSL_IMAG (gsldata) == 0.0)
		{
			GSL_SET_COMPLEX (&z,  PI_2, 0);
		}
		else
		{
			z = arctan().inverse().gsldata;
		}
		
		return CComplex(z);
	}

	// 双曲线函数
	CComplex CComplex::sinh() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;

		GSL_SET_COMPLEX (&z, ::sinh (R) * ::cos (I), ::cosh (R) * ::sin (I));

		return CComplex(z);
	}

	CComplex CComplex::cosh() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;

		GSL_SET_COMPLEX (&z, ::cosh (R) * ::cos (I), ::sinh (R) * ::sin (I));

		return CComplex(z);
	}

	CComplex CComplex::tanh() const
	{
		double R = GSL_REAL (gsldata), I = GSL_IMAG (gsldata);
		
		gsl_complex z;
		
		if (::fabs(R) < 1.0) 
		{
			double D = ::pow (::cos (I), 2.0) + ::pow (::sinh (R), 2.0);
			
			GSL_SET_COMPLEX (&z, ::sinh (R) * ::cosh (R) / D, 0.5 * ::sin (2 * I) / D);
		}
		else
		{
			double D = ::pow (::cos (I), 2.0) + ::pow (::sinh (R), 2.0);
			double F = 1 + ::pow (::cos (I) / ::sinh (R), 2.0);
			
			GSL_SET_COMPLEX (&z, 1.0 / (::tanh (R) * F), 0.5 * ::sin (2 * I) / D);
		}

		return CComplex(z);
	}
	CComplex CComplex::sech() const
	{
		return cosh().inverse();
	}

	CComplex CComplex::csch() const
	{
		return sinh().inverse();
	}

	CComplex CComplex::coth() const
	{
		return tanh().inverse();
	}

	// 反双曲线函数
	CComplex CComplex::arcsinh() const
	{
		
		return mul_imag(1.0).arcsin().mul_imag(-1.0);
	}

	CComplex CComplex::arccosh() const
	{
		return arccos().mul_imag(GSL_IMAG(gsldata) > 0 ? -1.0 : 1.0);
	}

	CComplex CComplex::arccosh(double x)
	{
		gsl_complex z;
		
		if (x >= 1)
		{
			GSL_SET_COMPLEX (&z,  CMath::acosh (x), 0);
		}
		else
		{
			if (x >= -1.0)
			{
				GSL_SET_COMPLEX (&z, 0,  acos (x));
			}
			else
			{
				GSL_SET_COMPLEX (&z,  CMath::acosh (-x),  PI);
			}
		}

		return CComplex(z);
	}

	CComplex CComplex::arctanh() const
	{
		if (GSL_IMAG (gsldata) == 0.0)
		{
			return arctanh(GSL_REAL (gsldata));
		}
		else
		{
			return mul_imag(1.0).arctan().mul_imag(-1.0);		
		}
	}

	CComplex CComplex::arctanh(double x)
	{
		gsl_complex z;
		
		if (x > -1.0 && x < 1.0)
		{
			GSL_SET_COMPLEX (&z, CMath::atanh (x), 0);
		}
		else
		{
			GSL_SET_COMPLEX (&z, CMath::atanh (1 / x), (x < 0) ?  PI_2 : -1* PI_2);
		}
		return CComplex(z);
	}

	CComplex CComplex::arcsech() const
	{
		return inverse().arccosh();
	}

	CComplex CComplex::arccsch() const
	{
		return inverse().arcsin();
	}

	CComplex CComplex::arccoth() const
	{
		return inverse().arctanh();
	}

	
}
