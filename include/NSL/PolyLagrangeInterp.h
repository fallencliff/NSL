/********************************************************************
	filename: 	PolyLagrangeInterp.h
	author:		hu zhijian
	created:	5:5:2010   16:53
	brief:	��Ԫ�������ղ�ֵ��
*********************************************************************/
#ifndef NSL_POLYLAGRANGEINTERP_H__
#define NSL_POLYLAGRANGEINTERP_H__

#include <NSL.h>
#include <BaseInterp.h>
#include <interp_struct.h>
#include <vector>


namespace gslcpp
{
	class NSL_EXPORT CPolyLagrangeInterp : public CBaseInterp  
	{
	public:

		CPolyLagrangeInterp(const interp_table& interp_data);

		virtual ~CPolyLagrangeInterp();

	private:	
		double RawInterp(const double* x, size_t jlo[]);

		double Mul(const double* x, double t, size_t m, size_t i); //����

		double Add(size_t next); //����

		std::vector<size_t> global_array;	//ÿ�����ӵ�ǰ���±�

		const double* x_to_interp; //������

		std::vector<size_t> index_base; //����ȷ������ֵ��λ��

	};
}

#endif // NSL_POLYLAGRANGEINTERP_H__
