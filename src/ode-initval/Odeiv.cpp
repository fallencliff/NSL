/********************************************************************
	filename: 	Odeiv.cpp
	author:		hu zhijian
	created:	7:5:2010   11:21
	brief:	solving ordinary differential equation (ODE) initial value problems
*********************************************************************/

#include <Odeiv.h>
#include <odeiv_type.h>
#include <ode-initval/ode_include.h>
namespace gslcpp
{
	const  COdeiv::StepTypeArray COdeiv::step_type_stack = COdeiv::StepTypeInit(); //初始化生成器的类型堆栈
	const  COdeiv::ControlTypeArray COdeiv::control_type_stack = COdeiv::ControlTypeInit(); //初始化生成器的类型堆栈

	COdeiv::COdeiv(const gsl_odeiv_system * ode_system,const Enum_Step_Type stepType /*= STEP_RK4*/,const Enum_Control_Type controlType /*= STANDARD*/ )
		:dydt(ode_system), s(NULL), c(NULL), e(NULL)
	{
		if (stepType > step_type_stack.size() || stepType < Enum_Step_Type(0))
		{
			gsl_error("type of step not exists",
				__FILE__, __LINE__, GSL_ENOMEM);
		}
		
		if (controlType > control_type_stack.size() || controlType < Enum_Control_Type(0))
		{
			gsl_error("type of control not exists",
				__FILE__, __LINE__, GSL_ENOMEM);
		}

		step_type = step_type_stack[stepType];
		control_type = control_type_stack[controlType];

		s = StepAlloc(ode_system->dimension);
		//c = ControlAlloc();

		e = Alloc(ode_system->dimension);
	}
	
	COdeiv::~COdeiv()
	{
		EvolveFree();
		ControlFree();
		StepFree();
	}
	
	COdeiv::StepTypeArray COdeiv::StepTypeInit()
	{
		COdeiv::StepTypeArray res;
 		res.push_back(gsl_odeiv_step_rk2);
 		res.push_back(gsl_odeiv_step_rk4);
 		res.push_back(gsl_odeiv_step_rkf45);
 		res.push_back(gsl_odeiv_step_rkck);
 		res.push_back(gsl_odeiv_step_rk8pd);
 		res.push_back(gsl_odeiv_step_rk2imp);
 		res.push_back(gsl_odeiv_step_rk4imp);
//		res.push_back(gsl_odeiv_step_bsimp);
 		res.push_back(gsl_odeiv_step_gear1);
// 		res.push_back(gsl_odeiv_step_gear2);
		return res;
	}
	
	COdeiv::ControlTypeArray COdeiv::ControlTypeInit()
	{
		COdeiv::ControlTypeArray res;

		res.push_back(gsl_odeiv_control_standard);
		res.push_back(gsl_odeiv_control_scaled);

		return res;
	}


	gsl_odeiv_step* COdeiv::StepAlloc( size_t dim )
	{
		gsl_odeiv_step *s = new  gsl_odeiv_step;
		
		if (s == 0)
		{
			GSL_ERROR_NULL ("failed to allocate space for ode struct", GSL_ENOMEM);
		};
		
		s->type = this->step_type;
		s->dimension = dim;
		
		s->state = s->type->alloc(dim);
		
		if (s->state == 0)
		{
			delete s;	 /* exception in constructor, avoid memory leak */
			
			s = NULL;
			
			GSL_ERROR_NULL ("failed to allocate space for ode state", GSL_ENOMEM);
		};
		
		return s;		
	}
	
	const char* COdeiv::StepName() const
	{
		return s->type->name;		
	}
	
	unsigned int COdeiv::StepOrder() const
	{
		return s->type->order(s->state);		
	}
	
	int COdeiv::StepReset()
	{
		return s->type->reset(s->state, s->dimension);	
	}
	
	void COdeiv::StepFree()
	{
		if (s)
		{
			s->type->free(s->state);
			delete s;
			s = NULL;
		}
	
	}
	
	int COdeiv::StepApply( double t, double h, double y[], double yerr[], const double dydt_in[], double dydt_out[], const gsl_odeiv_system * dydt )
	{
		return s->type->apply(s->state, s->dimension, t, h, y, yerr, dydt_in, dydt_out, dydt);		
	}
	
	gsl_odeiv_control * COdeiv::ControlAlloc()
	{
		gsl_odeiv_control * c =  new gsl_odeiv_control;
		
		if(c == 0) 
		{
			GSL_ERROR_NULL ("failed to allocate space for control struct", 
				GSL_ENOMEM);
		}
		
		c->type = control_type;
		c->state = c->type->alloc();
		
		if (c->state == 0)
		{
			delete c;       /* exception in constructor, avoid memory leak */
			c = NULL;
			GSL_ERROR_NULL ("failed to allocate space for control state", 
				GSL_ENOMEM);
		}
		
		return c;		
	}
	
	int COdeiv::ControlInit( double eps_abs, double eps_rel, double a_y, double a_dydt )
	{
		return c->type->init (c->state, eps_abs, eps_rel, a_y, a_dydt);		
	}
	
	void COdeiv::ControlFree()
	{
		if (c)
		{
			c->type->free(c->state);
			delete c;
			c = NULL;
		}
		
	}
	
	const char * COdeiv::ControlName() const
	{
		return c->type->name;		
	}
	
	int COdeiv::ControlHadjust( const double y0[], const double yerr[], const double dydt[], double * h )
	{
		return c->type->hadjust(c->state, s->dimension, s->type->order(s->state),
                          y0, yerr, dydt, h);		
	}
	
	int COdeiv::Control_Standard_New( double eps_abs, double eps_rel, double a_y, double a_dydt )
	{
		c = ControlAlloc ();
		
		int status = ControlInit (eps_abs, eps_rel, a_y, a_dydt);
		
		if (status != GSL_SUCCESS)
		{
			ControlFree();
			GSL_ERROR_NULL ("error trying to initialize control", status);
		}
		
		return GSL_SUCCESS;
	}
	
	int COdeiv::Reset()
	{
		e->count = 0;
		e->failed_steps = 0;
		e->last_step = 0.0;
		return GSL_SUCCESS;
	}
	
	void COdeiv::EvolveFree()
	{
		
		if (e)
		{
			delete [] e->dydt_out;
			delete [] e->dydt_in;
			delete [] e->yerr;
			delete [] e->y0;
			delete e;
			e = NULL;
		}
		
	}
	
	gsl_odeiv_evolve * COdeiv::Alloc( size_t dim )
	{
		gsl_odeiv_evolve *e = new gsl_odeiv_evolve;
		
		if (e == 0)
		{
			GSL_ERROR_NULL ("failed to allocate space for evolve struct",
				GSL_ENOMEM);
		}
		
		e->y0 = (double *) new double[dim];
		
		if (e->y0 == 0)
		{
			delete e;
			e = NULL;
			GSL_ERROR_NULL ("failed to allocate space for y0", GSL_ENOMEM);
		}
		
		e->yerr = (double *) new double[dim];
		
		if (e->yerr == 0)
		{
			delete [] e->y0;
			delete e;
			e = NULL;
			GSL_ERROR_NULL ("failed to allocate space for yerr", GSL_ENOMEM);
		}
		
		e->dydt_in = (double *) new double[dim];
		
		if (e->dydt_in == 0)
		{
			delete [] e->yerr;
			delete [] e->y0;
			delete e;
			e = NULL;
			GSL_ERROR_NULL ("failed to allocate space for dydt_in", GSL_ENOMEM);
		}
		
		e->dydt_out = (double *) new double[dim];
		
		if (e->dydt_out == 0)
		{
			delete [] e->dydt_in;
			delete [] e->yerr;
			delete [] e->y0;
			delete e;
			e = NULL;
			GSL_ERROR_NULL ("failed to allocate space for dydt_out", GSL_ENOMEM);
		}
		
		e->dimension = dim;
		e->count = 0;
		e->failed_steps = 0;
		e->last_step = 0.0;
		
		return e;		
	}
	
	int COdeiv::Apply( double&t, double t1, double& h, double y[] )
	{
		const double t0 = t;
		double h0 = h;
		int step_status;
		int final_step = 0;
		double dt = t1 - t0;  /* remaining time, possibly less than h */
		
		if (e->dimension != s->dimension)
		{
			GSL_ERROR ("step dimension must match evolution size", GSL_EINVAL);
		}
		
		if ((dt < 0.0 && h0 > 0.0) || (dt > 0.0 && h0 < 0.0))
		{
			GSL_ERROR ("step direction must match interval direction", GSL_EINVAL);
		}
		
		/* No need to copy if we cannot control the step size. */
		
		if (c != NULL)
		{
			DBL_MEMCPY (e->y0, y, e->dimension);
		}
		
		/* Calculate initial dydt once if the method can benefit. */
		
		if (s->type->can_use_dydt_in)
		{
			int status = GSL_ODEIV_FN_EVAL (dydt, t0, y, e->dydt_in);
			
			if (status) 
			{
				return status;
			}
		}
		
try_step:
		
		if ((dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt))
		{
			h0 = dt;
			final_step = 1;
		}
		else
		{
			final_step = 0;
		}
		
		if (s->type->can_use_dydt_in)
		{
			step_status = StepApply (t0, h0, y, e->yerr, e->dydt_in, e->dydt_out, dydt);
		}
		else
		{
			step_status = StepApply (t0, h0, y, e->yerr, NULL, e->dydt_out, dydt);
		}
		
		/* Check for stepper internal failure */
		
		if (step_status != GSL_SUCCESS) 
		{
			return step_status;
		}
		
		e->count++;
		e->last_step = h0;
		
		if (final_step)
		{
			t = t1;
		}
		else
		{
			t = t0 + h0;
		}
		
		if (c != NULL)
		{
			/* Check error and attempt to adjust the step. */
			const int hadjust_status  = ControlHadjust (y, e->yerr, e->dydt_out, &h0);
			
			if (hadjust_status == GSL_ODEIV_HADJ_DEC)
			{
				/* Step was decreased. Undo and go back to try again. */
				DBL_MEMCPY (y, e->y0, dydt->dimension);
				e->failed_steps++;
				goto try_step;
			}
		}
		
		h = h0;  /* suggest step size for next time-step */

		return step_status;		
	}
		
}