#ifndef  GSLCPP_CGCOMPLEX_H
#define  GSLCPP_CGCOMPLEX_H

/** 
 * @file Complex.h
 *
 * @brief ������
 *
 * @author 
 *
 * @date 2010-3-21
 */

#include <SL_dll.h>


#define GSL_REAL(z)     ((z).dat[0])
#define GSL_IMAG(z)     ((z).dat[1])
#define GSL_COMPLEX_P(zp) ((zp)->dat)
#define GSL_COMPLEX_P_REAL(zp)  ((zp)->dat[0])
#define GSL_COMPLEX_P_IMAG(zp)  ((zp)->dat[1])
#define GSL_COMPLEX_EQ(z1,z2) (((z1).dat[0] == (z2).dat[0]) && ((z1).dat[1] == (z2).dat[1]))

#define GSL_SET_COMPLEX(zp,x,y) do {(zp)->dat[0]=(x); (zp)->dat[1]=(y);} while(0)
#define GSL_SET_REAL(zp,x) do {(zp)->dat[0]=(x);} while(0)
#define GSL_SET_IMAG(zp,y) do {(zp)->dat[1]=(y);} while(0)

namespace gslcpp 
{

	/** 
	 * @brief ������
	 *
	 */
	class SL_DLL_API CComplex 
	{
		typedef struct
		{
			double dat[2];
		}
		gsl_complex;

		private:

			gsl_complex gsldata;	

		public:
			/** 
			 * @brief Ĭ�Ϲ��캯��, ���츴�� $z=0 + 0i$.
			 */
			CComplex();

			/** 
			 * @brief ���츴�� $z=r + i$.
			 */
			CComplex(const double r, const double i);

			/** 
			 * @brief �������캯��.
			 */
			CComplex(const CComplex& c);

			/** 
			 * @brief ���ι��캯��.
			 */
			CComplex(gsl_complex z);

			/** 
			 * @brief ���صѿ���ֱ������ϵ�µĸ���.	
			 * 
			 * @param ����xΪʵ��
			 * @param ����yΪ�鲿
			 * 
			 * @return ���ظ��� $z=x + yi$
			 *
			 */
			static CComplex rect(double x, double y);

			/** 
			 * @brief ���ؼ������µĸ���.	
			 * 
			 * @param r �뾶
			 * @param theta �Ƕ�
			 * 
			 * @return ���ظ��� $z=re^{i*theta} = r*(cos(theta)+i*sin(theta))$
			 */
			static CComplex polar(double r, double theta);

			/** 
			 * @brief ����ʵ��. ** @org: GSL_REAL
			 */
			double real() const;

			/** 
			 * @brief �����鲿.	
			 *
			 */
			double imag() const;

			/** 
			 * @brief �õѿ����������ø���	
			 * 
			 * @param x ʵ��
			 * @param y �鲿
			 *
			 */
			void set(double x, double y);

			/** 
			 * @brief ����ʵ��	
			 *
			 */
			void set_real(double a);

			/** 
			 * @brief �����鲿	
			 *
			 */
			void set_imag(double a);

			/** 
			 * @brief ���ط���	
			 *
			 */
			double arg() const;

			/** 
			 * @brief ����ģֵ. $|z|$		
			 *
			 */
			double abs() const;

			/** 
			 * @brief ����ģֵ��ƽ��. $|z|^2$		
			 *
			 */
			double abs2() const;

			/** 
			 * @brief ���� $log|z|$		
			 *
			 */
			double logabs() const;

			/*
			 * @brief ���ظ���z����10��logֵ, $log_10(z)$.	
			 *
			 */
			CComplex log10() const;

			/*
			 * @brief ���ظ���z����e��logֵ, $log_e(z)$.	
			 */
			CComplex log() const;

			/*
			 * @brief ���ظ���z���ڸ���b��logֵ, $log_b(z)$.	
			 */
			CComplex log_b(const CComplex& b) const;

			/*
			 * @brief $exp(z)$.	
			 */
			CComplex exp() const;
			/** 
			 * @brief ���ڸ���c,����һ������z_n, $z_n = z+c$	
			 *
			 */
			CComplex operator+(const CComplex& c) const;

			/** 
			 * @brief ���ڸ���c,����һ������z_n, $z_n = z-c$	
			 *
			 */
			CComplex operator-(const CComplex& c) const;	

			/** 
			 * @brief ���ڸ���c,����һ������z_n, $z_n = z*c$	
			 *
			 */
			CComplex operator*(const CComplex& c) const;

			/** 
			 * @brief ���ڸ���c,����һ������z_n, $z_n = z/c$	
			 *
			 */
			CComplex operator/(const CComplex& c) const;

			/** 
			 * @brief ����ʵ��x,����һ������z_n, $z_n = z+x$	
			 *
			 */
			CComplex operator+(double x) const;

			/** 
			 * @brief ����ʵ��x,����һ������z_n, $z_n = z-x$	
			 *
			 */
			CComplex operator-(double x) const;

			/** 
			 * @brief ����ʵ��x,����һ������z_n, $z_n = z*x$	
			 */
			CComplex operator*(double x) const;

			/** 
			 * @brief ����ʵ��x,����һ������z_n, $z_n = z/x$	
			 *
			 */
			CComplex operator/(double x) const;

			/** 
			 * @brief ����ʵ��x,����һ������z_n, $z_n = x+0*i$
			 *
			 */
			CComplex operator=(double x);

			/** 
			 * @brief ���ڸ���c,����һ������z_n, $z_n = z + c$
			 *
			 */
			CComplex operator=(const CComplex& c);

			/** 
			 * @brief �����븴����ȱȽ�
			 *
			 */
			int operator==(const CComplex& c) const;

			/** 
			 * @brief ����ʵ����ȱȽ�
			 *
			 */
			int CComplex::operator==(double x) const;

			/** 
			 * @brief ����ʵ��y,����һ������z_n, $z_n = z+iy$	
			 *
			 */
			CComplex add_imag(double y) const;

			/** 
			 * @brief ����ʵ��y,����һ������z_n, $z_n = z-iy$	
			 *
			 */
			CComplex sub_imag(double y) const;

			/** 
			 * @brief ����ʵ��y,����һ������z_n, $z_n = z*iy$	
			 *
			 */
			CComplex mul_imag(double y) const;

			/** 
			 * @brief ����ʵ��y,����һ������z_n, $z_n = z/iy$	
			 *
			 */
			CComplex div_imag(double y) const;

			/** 
			 * @brief ����z�Ĺ����, $z_n =x - iy$	
			 *
			 */
			CComplex conjugate() const;

			/** 
			 * @brief ����һ������z_n, $z_n = 1/z$	
			 *
			 */
			CComplex inverse() const;

			/** 
			 * @brief ����һ������z_n, $z_n = -z$	
			 *
			 */
			CComplex operator-(void) const;


			// ���Ǻ���

			/** 
			 * @brief ���ظ���z��sinֵ, $sin(z) = (e^{iz} - e^{-iz})/(2i)$.	
			 *
			 */
			CComplex sin() const;

			/** 
			 * @brief ���ظ���z��cosֵ, $cos(z) = (e^{iz} + e^{-iz})/2$.	
			 *
			 */
			CComplex cos() const;

			/** 
			 * @brief ���ظ���z��tanֵ, $tan(z) = \sin(z)/\cos(z)$.		
			 *
			 */
			CComplex tan() const;

			/** 
			 * @brief ���ظ���z��secֵ, $sec(z) = 1/\cos(z)$.	
			 *
			 */
			CComplex sec() const;

			/** 
			 * @brief ���ظ���z��cscֵ, $csc(z) = 1/\sin(z)$.	
			 *
			 */
			CComplex csc() const;

			/** 
			 * @brief ���ظ���z��cotֵ, $cot(z) = 1/\tan(z)$.  
			 *
			 */
			CComplex cot() const;


			// �����Ǻ���
			
			/** 
			 * @brief ���ظ���z��arcsinֵ, $arcsin(z)$.  
			 *
			 */
			CComplex arcsin() const;

			/** 
			 * @brief ����ʵ��x��arcsinֵ, $arcsin(x)$.  
			 *
			 */
			static CComplex arcsin(double x);

			/** 
			 * @brief ���ظ���z��arccosֵ, $arccos(z)$. 
			 *
			 */
			CComplex arccos() const;

			/** 
			 * @brief ����ʵ��x��arccosֵ, $arccos(x)$. 
			 *
			 */
			static CComplex arccos(double x);

			/** 
			 * @brief ���ظ���z��arctanֵ, $arctan(z)$. 
			 *
			 */
			CComplex arctan() const;

			/** 
			 * @brief ���ظ���z��arcsecֵ, $arcsec(z) = arccos(1/z)$. 
			 *
			 */
			CComplex arcsec() const;

			/** 
			 * @brief ����ʵ��x��arcsecֵ, $arcsec(x) = arccos(1/x)$. 
			 *
			 */
			static CComplex arcsec(double x);

			/** 
			 * @brief ���ظ���z��arccscֵ, $arccsc(z) = arcsin(1/z)$. 
			 *
			 */
			CComplex arccsc() const;

			/** 
			 * @brief ����ʵ��x��arccscֵ, $arccsc(x) = arcsin(1/x)$. 
			 *
			 */
			static CComplex arccsc(double x);

			/** 
			 * @brief ���ظ���z��arccotֵ, $arccot(z) = arctan(1/z)$. 
			 *
			 */
			CComplex arccot() const;

			// ˫���ߺ���

			/** 
			 * @brief ���ظ���z��sinhֵ, $sin(z) = (e^{z} - e^{-z})/2$. 
			 *
			 */
			CComplex sinh() const;

			/** 
			 * @brief ���ظ���z��coshֵ, $cos(z) = (e^{z} + e^{-z})/2$. 
			 *
			 */
			CComplex cosh() const;

			/** 
			 * @brief ���ظ���z��tanhֵ, $tanh(z) = sinh(z)/cosh(z)$. 
			 *
			 */
			CComplex tanh() const;

			/** 
			 * @brief ���ظ���z��sechֵ, $sech(z) = 1/cosh(z)$. 
			 *
			 */
			CComplex sech() const;

			/** 
			 * @brief ���ظ���z��cschֵ, $csch(z) = 1/sinh(z)$.	
			 *
			 */
			CComplex csch() const;

			/** 
			 * @brief ���ظ���z��cothֵ, $coth(z) = 1/tanh(z)$.	
			 *
			 */
			CComplex coth() const;

			// ��˫���ߺ���
			
			/** 
			 * @brief ���ظ���z��arcsinhֵ, $arcsinh(z)$.	
			 *
			 */
			CComplex arcsinh() const;

			/** 
			 * @brief ���ظ���z��arccoshֵ, $arccosh(z)$.	
			 *
			 */
			CComplex arccosh() const;

			/** 
			 * @brief ����ʵ��x��arccoshֵ, $arccosh(x)$.	
			 *
			 */
			static CComplex arccosh(double x);

			/** 
			 * @brief ���ظ���z��arctanhֵ, $arctanh(z)$.	
			 *
			 */
			CComplex arctanh() const;

			/** 
			 * @brief ����ʵ��x��arctanhֵ, $arctanh(x)$.	
			 *
			 */
			static CComplex arctanh(double x);

			/** 
			 * @brief ���ظ���z��arcsechֵ, $arcsech(z) = arccosh(1/z)$.	
			 *
			 */
			CComplex arcsech() const;

			/** 
			 * @brief ���ظ���z��arccschֵ, $arccsch(z) = arcsinh(1/z)$.	
			 *
			 */
			CComplex arccsch() const;
			
			/** 
			 * @brief ���ظ���z��arccothֵ, $arccoth(z) = arctanh(1/z)$.	
			 *
			 */
			CComplex arccoth() const;
			
			/*
			 * @brief ���ڸ���c, ����һ������z_n, $z_n = z^c$.	
			 *
			 */
			CComplex pow(CComplex& c) const;
			
			/*
			 * @brief ����ʵ��x, ����һ������z_n, $z_n = z^x$.	
			 *
			 */
			CComplex pow(double x) const;

			/*
			 * @brief ���ظ���zƽ����,ʵ�����ڵ�����, $sqrt(z)$.	
			 *
			 */
			CComplex sqrt();
			
			/*
			 * @brief ����ʵ��xƽ����,ʵ�����ڵ����� $sqrt(x)$.	
			 *
			 */
			static CComplex sqrt(double x);
	};
}

#endif
