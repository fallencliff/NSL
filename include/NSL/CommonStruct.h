/********************************************************************
	filename: 	CommonStruct.h
	author:		hu zhijian
	created:	13:5:2010   17:59
	brief:	
*********************************************************************/
#ifndef NSL_COMMONSTRUCT_H__
#define NSL_COMMONSTRUCT_H__
#include <stdio.h>

namespace gslcpp
{
	class CBlock;

	typedef struct
	{
		size_t size;
		double *data;
	}gsl_block;

	typedef struct 
	{
		size_t size;
		size_t stride;
		double *data;
		CBlock *block;
		int owner;
	}gsl_vector;

	typedef struct 
	{
		size_t size1;
		size_t size2;
		size_t tda;
		double * data;
		CBlock * block;
		int owner;
	}gsl_matrix;

	typedef struct
	{
		size_t n;
		size_t k;
		size_t *data;
	}gsl_combination;

	typedef struct
	{
		double dat[2];
	}gsl_complex;


	/************************************************************************/
	/* for CHistogram
	/************************************************************************/
	typedef struct 
	{
		size_t n ;
		double * range ;
		double * bin ;
	} gsl_histogram;
	
	typedef struct 
	{
		size_t n ;
		double * range ;
		double * sum ;
	}gsl_histogram_pdf;

	/************************************************************************/
	/* for CNtuple  
	/************************************************************************/
	typedef struct 
	{
		FILE * file;
		void * ntuple_data;
		size_t size;
		
	} gsl_ntuple;
	
	typedef struct 
	{
		int (* function) (void * ntuple_data, void * params);
		void * params;
	}select_fn;
	
	typedef struct 
	{
		double (* function) (void * ntuple_data, void * params);
		void * params;
	}value_fn;


	/************************************************************************/
	/* for CPermutation
	/************************************************************************/
	typedef struct
	{
		size_t size;
		size_t *data;
	}gsl_permutation;
}






#endif // NSL_COMMONSTRUCT_H__