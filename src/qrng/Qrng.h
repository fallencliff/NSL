/********************************************************************
	filename: 	Qrng.h
	author:		hu zhijian
	created:	5:5:2010   16:42
	brief:	准随机数生成器
*********************************************************************/

#ifndef NWPU_QRNG_H
#define NWPU_QRNG_H
#pragma warning( disable : 4251) //type_stack有警告

#include <SL_DLL.h>
#include <Qrng_type.h>
#include <vector>
using std::vector;
namespace gslcpp
{
	class SL_DLL_API CQrng
	{
	public:
		typedef struct
		{
			const char * name;
			unsigned int max_dimension;
			size_t (*state_size) (unsigned int dimension);
			int (*init_state) (void * state, unsigned int dimension);
			int (*get) (void * state, unsigned int dimension, double x[]);
		}gsl_qrng_type;

		typedef vector<const gsl_qrng_type*> TypeArray;
		
		static const TypeArray type_stack;
		
		static CQrng::TypeArray TypeInit();	

	private:
		const gsl_qrng_type * qrng_type;
		unsigned int dimension;
		size_t state_size;
		void * state;


	public:
		CQrng(const size_t dimension, QRNG_ORDER type = SOBOL);

		virtual ~CQrng();

		const char* Name() const;

		size_t StateSize() const;

		void* State() const;

		int Get(double x[], size_t size) const;

		size_t MaxDimension() const;

		const gsl_qrng_type* GetType() const;

		CQrng& operator=(const CQrng& other);
	};



}



#endif