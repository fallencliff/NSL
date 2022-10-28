/********************************************************************
	filename: 	odeiv_struct.h
	author:		hu zhijian
	created:	7:5:2010   11:34
	brief:	struct declaration for COdeiv
*********************************************************************/
#ifndef NWPU_ODEIV_STRUCT_H__
#define NWPU_ODEIV_STRUCT_H__

namespace gslcpp
{
#define GSL_ODEIV_FN_EVAL(S,t,y,f)  (*((S)->function))(t,y,f,(S)->params)
#define GSL_ODEIV_JA_EVAL(S,t,y,dfdy,dfdt)  (*((S)->jacobian))(t,y,dfdy,dfdt,(S)->params)

/* Possible return values for an hadjust() evolution method.
	*/
#define GSL_ODEIV_HADJ_INC   1  /* step was increased */
#define GSL_ODEIV_HADJ_NIL   0  /* step unchanged     */
#define GSL_ODEIV_HADJ_DEC (-1) /* step decreased     */
	
	typedef struct  
	{
		int (* function) (double t, const double y[], double dydt[], void * params);
		int (* jacobian) (double t, const double y[], double * dfdy, double dfdt[], void * params);
		size_t dimension;
		void * params;
	}gsl_odeiv_system;

	typedef struct 
	{
		const char * name;
		int can_use_dydt_in;
		int gives_exact_dydt_out;
		void * (*alloc) (size_t dim);
		int  (*apply)  (void * state, size_t dim, double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt);
		int  (*reset) (void * state, size_t dim);
		unsigned int  (*order) (void * state);
		void (*free)  (void * state);
	}gsl_odeiv_step_type;
	
	typedef struct 
	{
		const gsl_odeiv_step_type * type;
		size_t dimension;
		void * state;
	}gsl_odeiv_step;

	typedef struct
	{
		const char * name;
		void * (*alloc) (void);
		int  (*init) (void * state, double eps_abs, double eps_rel, double a_y, double a_dydt);
		int  (*hadjust) (void * state, size_t dim, unsigned int ord, const double y[], const double yerr[], const double yp[], double * h);
		void (*free) (void * state);
	}gsl_odeiv_control_type;
	
	typedef struct
	{
		const gsl_odeiv_control_type * type;
		void * state;
	}gsl_odeiv_control;

	typedef struct 
	{
		size_t dimension;
		double * y0;
		double * yerr;
		double * dydt_in;
		double * dydt_out;
		double last_step;
		unsigned long int count;
		unsigned long int failed_steps;
	}gsl_odeiv_evolve;
}





#endif // NWPU_ODEIV_STRUCT_H__

