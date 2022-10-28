/********************************************************************
	filename: 	Qrng.cpp
	author:		huzhijian
	created:	5:5:2010   16:42
	brief:	
*********************************************************************/
#include <gsl_errno.h>
#include <Qrng.h>
#include <qrng/QRandAlgorithm.h>


namespace gslcpp
{
	const  CQrng::TypeArray CQrng::type_stack = CQrng::TypeInit(); //初始化生成器的类型堆栈
	
	CQrng::TypeArray CQrng::TypeInit()
	{
		CQrng::TypeArray res;
		res.push_back(&nied2_type);
		res.push_back(&sobol_type);
		return res;
	}
	CQrng::CQrng(const size_t dimension, QRNG_ORDER type /*= SOBOL*/)
		:state(NULL), dimension(0), qrng_type(NULL)
	{
		if (type > type_stack.size() || type < QRNG_ORDER(0))
		{
			gsl_error("type of generator not exists",
				__FILE__, __LINE__, GSL_ENOMEM);
		}

		qrng_type = type_stack[type];

		this->dimension = dimension;

		this->state_size = qrng_type->state_size(dimension);

		/* FIXME: we're assuming that a char is 8 bits */
		this->state = new char[this->state_size];
	
		if (state == NULL)
		{		
			gsl_error("failed to allocate space for qrng state",
				__FILE__, __LINE__, GSL_ENOMEM);
		}
		
		qrng_type->init_state(this->state, dimension);
	
	}
	
	CQrng::~CQrng()
	{
		if (state)
		{
			delete [] (char*)state;
		}
			
	}
	
	inline const char* CQrng::Name() const
	{
		return qrng_type->name;		
	}
	
	inline size_t CQrng::StateSize() const
	{
		return state_size;		
	}
	
	inline void* CQrng::State() const
	{
		return state;	
	}
	
	inline int CQrng::Get( double x[], size_t size ) const
	{
		if (size != dimension)
		{
			gsl_error("The space available for x must match the dimension of the generator",
				__FILE__, __LINE__, GSL_EINVAL);
		}
		return (qrng_type->get) (state, dimension, x);		
	}
	
	inline size_t CQrng::MaxDimension() const
	{
		return 	qrng_type->max_dimension;	
	}
	
	CQrng& CQrng::operator=( const CQrng& other )
	{
		if (GetType() != other.GetType())
		{
			gsl_error ("generators must be of the same type",
                        __FILE__, __LINE__, GSL_EINVAL);
		}
		dimension = other.dimension;
		state_size = other.state_size;
		memcpy (state, other.state, other.state_size);

		return *this;
	}
	
	inline const gsl_qrng_type* CQrng::GetType() const
	{
		return qrng_type;		
	}
	

}