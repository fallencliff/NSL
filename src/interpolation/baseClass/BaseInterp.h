/********************************************************************
	filename: 	BaseInterp.h
	author:		hu zhijian
	created:	5:5:2010   16:46
	brief:	��ֵ��������,����������е����飬���˲�ֵ!
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
			bool cor; //�Ƿ�ʹ�ü�����
			
			size_t cache; //�ϴ�Ѱ�ҵ��±�
			
			size_t dj; // ������ֶ�
		}ACCEL;

	private:
		void Free();
		void Init(const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, const size_t& y_size);

	//protected:
		//std::vector<const CVector*> xx; //����Ա���
		
		std::vector<const double*> xx;//����Ա���
		const double* yy;//��������
		const size_t* xx_size;
		const size_t yy_size;
		const size_t dim; //ά��
		const size_t mm; //��ֵ�����ĸ���

		ACCEL* acc; //������
		size_t* jlo; //Ѱ�ҵ����±�,���ظ�������
		bool extroInterp_flag;

	public:	
		CBaseInterp():acc(NULL), jlo(NULL), xx_size(0), yy(NULL), yy_size(0), dim(0), mm(0)
		{

		}
		CBaseInterp(const size_t& dimension, const double* const* x_Value, const size_t* x_data_size, 
			const double* y_Value, const size_t& y_size, const size_t& m, bool extro_flag);
		
		virtual ~CBaseInterp();
		
		//�ڵ�iά�Ա�����Ѱ��x����Ӧ���±꣬ȫȡ�����ַ�
		size_t Locate(const double& x, const size_t& i);
		
		//�ڵ�iά�Ա�����Ѱ��x����Ӧ���±�, ʹ�ü����� ,�ֲ����ַ�
		size_t Hunt(const double& x, const size_t& i);
		
		//�û��ӿ�
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

		//�����������д�˺��������ز�ֵ���
		double virtual RawInterp(const double* x, size_t jlo[]) = 0;
	};
}

#endif // BASEINTERP_H__
