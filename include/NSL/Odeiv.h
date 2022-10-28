/********************************************************************
	filename: 	Odeiv.h
	author:		hu zhijian
	created:	7:5:2010   11:19
	brief:	solving ordinary differential equation (ODE) initial value problems
*********************************************************************/
#ifndef NSL_ODEIV_H__
#define NSL_ODEIV_H__


#pragma warning( disable : 4251) //type_stackÓÐ¾¯¸æ

#include <NSL.h>
#include <odeiv_type.h>
#include <odeiv_struct.h>
#include <Vector>

using std::vector;

namespace gslcpp
{
	class NSL_EXPORT COdeiv
	{
	public:

		typedef vector<const gsl_odeiv_step_type*> StepTypeArray;
		typedef vector<const gsl_odeiv_control_type*> ControlTypeArray;

		static const StepTypeArray step_type_stack;
		static const ControlTypeArray control_type_stack;

		static COdeiv::StepTypeArray StepTypeInit();
		static COdeiv::ControlTypeArray ControlTypeInit();

	private:
		//variable
		const gsl_odeiv_step_type *step_type;
		const gsl_odeiv_control_type * control_type;
		const gsl_odeiv_system * dydt;
		gsl_odeiv_step* s;
		gsl_odeiv_control* c;
		gsl_odeiv_evolve * e;

		//private functions
		gsl_odeiv_step* StepAlloc(size_t dim);
		gsl_odeiv_control * ControlAlloc();

		//step functions
		
		unsigned int StepOrder() const;
		int StepReset();
		void StepFree();
		int StepApply
			(
			double t,
			double h,
			double y[],
			double yerr[],
			const double dydt_in[],
			double dydt_out[],
			const gsl_odeiv_system * dydt
			);
		
		//control functions
		int ControlInit(double eps_abs, double eps_rel, double a_y, double a_dydt);
		void ControlFree();
		
		int ControlHadjust(const double y0[], const double yerr[], const double dydt[], double * h);

		//Evolution functions
		void EvolveFree();
	public:
		COdeiv(const gsl_odeiv_system * ode_system,
			   const Enum_Step_Type stepType = STEP_RK4,
			   const Enum_Control_Type controlType = STANDARD
			   );
		virtual ~COdeiv();

		int Reset();
		gsl_odeiv_evolve * Alloc (size_t dim);

		int Apply(double&t, double t1, double& h, double y[]);
		int Control_Standard_New(double eps_abs, double eps_rel, double a_y, double a_dydt);

		const char * ControlName() const;
		const char* StepName() const;
	};

}

#endif // NSL_ODEIV_H__