/********************************************************************
	filename: 	Qrng.h
	author:		hu zhijian
	created:	5:5:2010   16:42
	brief:	准随机数生成器
*********************************************************************/

#ifndef NSL_QRNG_H__
#define NSL_QRNG_H__

#pragma warning( disable : 4251) //type_stack dll interface 有警告

#include <NSL.h>
#include <Qrng_type.h>
#include <vector>

using std::vector;

namespace gslcpp
{
	class NSL_EXPORT CQrng
	{
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

#endif // NSL_QRNG_H__