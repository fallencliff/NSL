/********************************************************************
	filename: 	state.h
	author:		hu zhijian
	created:	6:5:2010   10:05
	brief:	miser, vegas, and plain state struct
*********************************************************************/
#ifndef NWPU_STATE_H__
#define NWPU_STATE_H__


#include <stdio.h>

namespace gslcpp
{
	typedef	struct
	{
		double (*f)(double * x_array, size_t dim, void * params);
		void* params;
	}MonteCarlo_function;

	typedef  struct
	{
		size_t dim;
		double *x;
	}monte_plain_state;

	typedef  struct 
	{
		size_t min_calls;
		size_t min_calls_per_bisection;
		double dither;
		double estimate_frac;
		double alpha;
		size_t dim;
		int estimate_style;
		int depth;
		int verbose;
		double * x;
		double * xmid;
		double * sigma_l;
		double * sigma_r;
		double * fmax_l;
		double * fmax_r;
		double * fmin_l;
		double * fmin_r;
		double * fsum_l;
		double * fsum_r;
		double * fsum2_l;
		double * fsum2_r;
		size_t * hits_l;
		size_t * hits_r;
	} monte_miser_state; 

	enum 
	{
		VEGAS_MODE_IMPORTANCE = 1, 
		VEGAS_MODE_IMPORTANCE_ONLY = 0, 
		VEGAS_MODE_STRATIFIED = -1
	};

	
	typedef  struct
	{
		/* grid */
		size_t dim;
		size_t bins_max;
		unsigned int bins;
		unsigned int boxes; /* these are both counted along the axes */
		double * xi;
		double * xin;
		double * delx;
		double * weight;
		double vol;
		
		double * x;
		int * bin;
		int * box;
		
		/* distribution */
		double * d;
		
		/* control variables */
		double alpha;
		int mode;
		int verbose;
		unsigned int iterations;
		int stage;
		
		/* scratch variables preserved between calls to vegas1/2/3  */
		double jac;
		double wtd_int_sum; 
		double sum_wgts;
		double chi_sum;
		double chisq;
		
		double result;
		double sigma;
		
		unsigned int it_start;
		unsigned int it_num;
		unsigned int samples;
		unsigned int calls_per_box;
		
		FILE * ostream;
		
	} monte_vegas_state;

}

#endif // NWPU_STATE_H__