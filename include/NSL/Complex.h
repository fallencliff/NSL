/********************************************************************
	filename: 	Complex.h
	author:		hu zhijian
	created:	5:5:2010   16:50
	brief:	复数类
*********************************************************************/

#ifndef NSL_COMPLEX_H__
#define NSL_COMPLEX_H__

#include <NSL.h>
#include <CommonStruct.h>

namespace gslcpp 
{

	/** 
	 * @brief 复数类
	 *
	 */
	class NSL_EXPORT CComplex 
	{
		private:

			gsl_complex gsldata;	

		public:
			/** 
			 * @brief 默认构造函数, 构造复数 $z=0 + 0i$.
			 */
			CComplex();

			/** 
			 * @brief 构造复数 $z=r + i$.
			 */
			CComplex(const double r, const double i);

			/** 
			 * @brief 拷贝构造函数.
			 */
			CComplex(const CComplex& c);

			/** 
			 * @brief 带参构造函数.
			 */
			CComplex(gsl_complex z);

			/** 
			 * @brief 返回笛卡尔直角坐标系下的复数.	
			 * 
			 * @param 参数x为实部
			 * @param 参数y为虚部
			 * 
			 * @return 返回复数 $z=x + yi$
			 *
			 */
			static CComplex rect(double x, double y);

			/** 
			 * @brief 返回极坐标下的复数.	
			 * 
			 * @param r 半径
			 * @param theta 角度
			 * 
			 * @return 返回复数 $z=re^{i*theta} = r*(cos(theta)+i*sin(theta))$
			 */
			static CComplex polar(double r, double theta);

			/** 
			 * @brief 返回实部. ** @org: GSL_REAL
			 */
			double real() const;

			/** 
			 * @brief 返回虚部.	
			 *
			 */
			double imag() const;

			/** 
			 * @brief 用笛卡尔坐标设置复数	
			 * 
			 * @param x 实部
			 * @param y 虚部
			 *
			 */
			void set(double x, double y);

			/** 
			 * @brief 设置实部	
			 *
			 */
			void set_real(double a);

			/** 
			 * @brief 设置虚部	
			 *
			 */
			void set_imag(double a);

			/** 
			 * @brief 返回幅角	
			 *
			 */
			double arg() const;

			/** 
			 * @brief 返回模值. $|z|$		
			 *
			 */
			double abs() const;

			/** 
			 * @brief 返回模值的平方. $|z|^2$		
			 *
			 */
			double abs2() const;

			/** 
			 * @brief 返回 $log|z|$		
			 *
			 */
			double logabs() const;

			/*
			 * @brief 返回复数z基于10的log值, $log_10(z)$.	
			 *
			 */
			CComplex log10() const;

			/*
			 * @brief 返回复数z基于e的log值, $log_e(z)$.	
			 */
			CComplex log() const;

			/*
			 * @brief 返回复数z基于复数b的log值, $log_b(z)$.	
			 */
			CComplex log_b(const CComplex& b) const;

			/*
			 * @brief $exp(z)$.	
			 */
			CComplex exp() const;
			/** 
			 * @brief 对于复数c,返回一个复数z_n, $z_n = z+c$	
			 *
			 */
			CComplex operator+(const CComplex& c) const;

			/** 
			 * @brief 对于复数c,返回一个复数z_n, $z_n = z-c$	
			 *
			 */
			CComplex operator-(const CComplex& c) const;	

			/** 
			 * @brief 对于复数c,返回一个复数z_n, $z_n = z*c$	
			 *
			 */
			CComplex operator*(const CComplex& c) const;

			/** 
			 * @brief 对于复数c,返回一个复数z_n, $z_n = z/c$	
			 *
			 */
			CComplex operator/(const CComplex& c) const;

			/** 
			 * @brief 对于实数x,返回一个复数z_n, $z_n = z+x$	
			 *
			 */
			CComplex operator+(double x) const;

			/** 
			 * @brief 对于实数x,返回一个复数z_n, $z_n = z-x$	
			 *
			 */
			CComplex operator-(double x) const;

			/** 
			 * @brief 对于实数x,返回一个复数z_n, $z_n = z*x$	
			 */
			CComplex operator*(double x) const;

			/** 
			 * @brief 对于实数x,返回一个复数z_n, $z_n = z/x$	
			 *
			 */
			CComplex operator/(double x) const;

			/** 
			 * @brief 对于实数x,返回一个复数z_n, $z_n = x+0*i$
			 *
			 */
			CComplex operator=(double x);

			/** 
			 * @brief 对于复数c,返回一个复数z_n, $z_n = z + c$
			 *
			 */
			CComplex operator=(const CComplex& c);

			/** 
			 * @brief 复数与复数相等比较
			 *
			 */
			int operator==(const CComplex& c) const;

			/** 
			 * @brief 复数实数相等比较
			 *
			 */
			int CComplex::operator==(double x) const;

			/** 
			 * @brief 对于实数y,返回一个复数z_n, $z_n = z+iy$	
			 *
			 */
			CComplex add_imag(double y) const;

			/** 
			 * @brief 对于实数y,返回一个复数z_n, $z_n = z-iy$	
			 *
			 */
			CComplex sub_imag(double y) const;

			/** 
			 * @brief 对于实数y,返回一个复数z_n, $z_n = z*iy$	
			 *
			 */
			CComplex mul_imag(double y) const;

			/** 
			 * @brief 对于实数y,返回一个复数z_n, $z_n = z/iy$	
			 *
			 */
			CComplex div_imag(double y) const;

			/** 
			 * @brief 返回z的共轭复数, $z_n =x - iy$	
			 *
			 */
			CComplex conjugate() const;

			/** 
			 * @brief 返回一个复数z_n, $z_n = 1/z$	
			 *
			 */
			CComplex inverse() const;

			/** 
			 * @brief 返回一个复数z_n, $z_n = -z$	
			 *
			 */
			CComplex operator-(void) const;


			// 三角函数

			/** 
			 * @brief 返回复数z的sin值, $sin(z) = (e^{iz} - e^{-iz})/(2i)$.	
			 *
			 */
			CComplex sin() const;

			/** 
			 * @brief 返回复数z的cos值, $cos(z) = (e^{iz} + e^{-iz})/2$.	
			 *
			 */
			CComplex cos() const;

			/** 
			 * @brief 返回复数z的tan值, $tan(z) = \sin(z)/\cos(z)$.		
			 *
			 */
			CComplex tan() const;

			/** 
			 * @brief 返回复数z的sec值, $sec(z) = 1/\cos(z)$.	
			 *
			 */
			CComplex sec() const;

			/** 
			 * @brief 返回复数z的csc值, $csc(z) = 1/\sin(z)$.	
			 *
			 */
			CComplex csc() const;

			/** 
			 * @brief 返回复数z的cot值, $cot(z) = 1/\tan(z)$.  
			 *
			 */
			CComplex cot() const;


			// 反三角函数
			
			/** 
			 * @brief 返回复数z的arcsin值, $arcsin(z)$.  
			 *
			 */
			CComplex arcsin() const;

			/** 
			 * @brief 返回实数x的arcsin值, $arcsin(x)$.  
			 *
			 */
			static CComplex arcsin(double x);

			/** 
			 * @brief 返回复数z的arccos值, $arccos(z)$. 
			 *
			 */
			CComplex arccos() const;

			/** 
			 * @brief 返回实数x的arccos值, $arccos(x)$. 
			 *
			 */
			static CComplex arccos(double x);

			/** 
			 * @brief 返回复数z的arctan值, $arctan(z)$. 
			 *
			 */
			CComplex arctan() const;

			/** 
			 * @brief 返回复数z的arcsec值, $arcsec(z) = arccos(1/z)$. 
			 *
			 */
			CComplex arcsec() const;

			/** 
			 * @brief 返回实数x的arcsec值, $arcsec(x) = arccos(1/x)$. 
			 *
			 */
			static CComplex arcsec(double x);

			/** 
			 * @brief 返回复数z的arccsc值, $arccsc(z) = arcsin(1/z)$. 
			 *
			 */
			CComplex arccsc() const;

			/** 
			 * @brief 返回实数x的arccsc值, $arccsc(x) = arcsin(1/x)$. 
			 *
			 */
			static CComplex arccsc(double x);

			/** 
			 * @brief 返回复数z的arccot值, $arccot(z) = arctan(1/z)$. 
			 *
			 */
			CComplex arccot() const;

			// 双曲线函数

			/** 
			 * @brief 返回复数z的sinh值, $sin(z) = (e^{z} - e^{-z})/2$. 
			 *
			 */
			CComplex sinh() const;

			/** 
			 * @brief 返回复数z的cosh值, $cos(z) = (e^{z} + e^{-z})/2$. 
			 *
			 */
			CComplex cosh() const;

			/** 
			 * @brief 返回复数z的tanh值, $tanh(z) = sinh(z)/cosh(z)$. 
			 *
			 */
			CComplex tanh() const;

			/** 
			 * @brief 返回复数z的sech值, $sech(z) = 1/cosh(z)$. 
			 *
			 */
			CComplex sech() const;

			/** 
			 * @brief 返回复数z的csch值, $csch(z) = 1/sinh(z)$.	
			 *
			 */
			CComplex csch() const;

			/** 
			 * @brief 返回复数z的coth值, $coth(z) = 1/tanh(z)$.	
			 *
			 */
			CComplex coth() const;

			// 反双曲线函数
			
			/** 
			 * @brief 返回复数z的arcsinh值, $arcsinh(z)$.	
			 *
			 */
			CComplex arcsinh() const;

			/** 
			 * @brief 返回复数z的arccosh值, $arccosh(z)$.	
			 *
			 */
			CComplex arccosh() const;

			/** 
			 * @brief 返回实数x的arccosh值, $arccosh(x)$.	
			 *
			 */
			static CComplex arccosh(double x);

			/** 
			 * @brief 返回复数z的arctanh值, $arctanh(z)$.	
			 *
			 */
			CComplex arctanh() const;

			/** 
			 * @brief 返回实数x的arctanh值, $arctanh(x)$.	
			 *
			 */
			static CComplex arctanh(double x);

			/** 
			 * @brief 返回复数z的arcsech值, $arcsech(z) = arccosh(1/z)$.	
			 *
			 */
			CComplex arcsech() const;

			/** 
			 * @brief 返回复数z的arccsch值, $arccsch(z) = arcsinh(1/z)$.	
			 *
			 */
			CComplex arccsch() const;
			
			/** 
			 * @brief 返回复数z的arccoth值, $arccoth(z) = arctanh(1/z)$.	
			 *
			 */
			CComplex arccoth() const;
			
			/*
			 * @brief 对于复数c, 返回一个复数z_n, $z_n = z^c$.	
			 *
			 */
			CComplex pow(CComplex& c) const;
			
			/*
			 * @brief 对于实数x, 返回一个复数z_n, $z_n = z^x$.	
			 *
			 */
			CComplex pow(double x) const;

			/*
			 * @brief 返回复数z平方根,实部大于等于零, $sqrt(z)$.	
			 *
			 */
			CComplex sqrt();
			
			/*
			 * @brief 返回实数x平方根,实部大于等于零 $sqrt(x)$.	
			 *
			 */
			static CComplex sqrt(double x);
	};
}

#endif // NSL_COMPLEX_H__
