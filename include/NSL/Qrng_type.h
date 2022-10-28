/********************************************************************
	filename: 	Qrng_type.h
	author:		hu zhijian
	created:	12:5:2010   19:43
	brief:	enumeration for different Quasi-Random algorithm
*********************************************************************/



#ifndef NSL_QRNG_TYPE_H__
#define NSL_QRNG_TYPE_H__


namespace gslcpp
{
	typedef enum
	{
		NIEDERREITER_BASE_2,
			SOBOL
	}QRNG_ORDER;
	
	//typedef enum QRNG_ORDER_T QRNG_ORDER;

	
	typedef struct
	{
		const char * name;
		unsigned int max_dimension;
		size_t (*state_size) (unsigned int dimension);
		int (*init_state) (void * state, unsigned int dimension);
		int (*get) (void * state, unsigned int dimension, double x[]);
		
	}gsl_qrng_type;
		
}

#endif // NSL_QRNG_TYPE_H__