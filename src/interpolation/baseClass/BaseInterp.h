/********************************************************************
	filename: 	BaseInterp.h
	author:		hu zhijian
	created:	5:5:2010   16:46
	brief:	插值基类声明,本类干了所有的事情，除了插值!
*********************************************************************/

#ifndef NWPU_BASEINTERP_H__
#define NWPU_BASEINTERP_H__


#pragma warning( disable : 4251)

#include <SL_dll.h>

#include <vector>

namespace gslcpp
{
	
	class SL_DLL_API CBaseInterp  
	{
	public:
		typedef struct
		{
			bool cor; //是否使用加速器
			
			size_t cache; //上次寻找的下标
			
			size_t dj; // 总区间分段
		}ACCEL;

	private:
		void Free();
		void Init(const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, const size_t& y_size);

	//protected:
		//std::vector<const CVector*> xx; //存放自变量
		
		std::vector<const double*> xx;//存放自变量
		const double* yy;//存放因变量
		const size_t* xx_size;
		const size_t yy_size;
		const size_t dim; //维数
		const size_t mm; //插值所需点的个数

		ACCEL* acc; //加速器
		size_t* jlo; //寻找到的下标,返回给派生类
		bool extroInterp_flag;

	public:	
		CBaseInterp():acc(NULL), jlo(NULL), xx_size(0), yy(NULL), yy_size(0), dim(0), mm(0)
		{

		}
		CBaseInterp(const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, 
			const double* y_Value, const size_t& y_size, const size_t& m, bool extro_flag);
		
		virtual ~CBaseInterp();
		
		//在第i维自变量中寻找x所对应的下标，全取决二分法
		size_t Locate(const double& x, const size_t& i);
		
		//在第i维自变量中寻找x所对应的下标, 使用加速器 ,局部二分法
		size_t Hunt(const double& x, const size_t& i);
		
		//用户接口
		double Interp(const double* x, const size_t& size);

		size_t GetNumOfDimension() const
		{
			return dim;
		}

		const double* Get_X_Addr(size_t i) const
		{
			return xx[i];
		}

		const double* Get_Y_Addr() const
		{
			return yy;
		}

		size_t Get_X_Size(const size_t& i) const 
		{
			return xx_size[i];
		}

		size_t Get_Y_Size() const 
		{
			return yy_size;
		}

		size_t Get_mm() const 
		{
			return mm;
		}

		const size_t* Get_jlo_Addr() const 
		{
			return jlo;
		}

		//派生类必须重写此函数，返回插值结果
		double virtual RawInterp(const double* x, size_t jlo[]) = 0;
	};
}

#endif // BASEINTERP_H__
