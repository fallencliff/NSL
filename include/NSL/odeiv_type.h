/********************************************************************
	filename: 	odeiv_type.h
	author:		hu zhijian
	created:	7:5:2010   19:42
	brief:	odeiv_step_type for COdeiv
*********************************************************************/

#ifndef NWPU_ODEIV_STEP_TYPE_H__
#define NWPU_ODEIV_STEP_TYPE_H__

typedef enum
{
	
	/************************************************************************/
	/* Embedded Runge-Kutta (2, 3) method. 
	/************************************************************************/
	STEP_RK2 = 0, 

	/************************************************************************/
	/* 4th order (classical) Runge-Kutta.
	/************************************************************************/
	STEP_RK4,

	/************************************************************************/
	/* Embedded Runge-Kutta-Fehlberg (4, 5) method. This method is a good generalpurpose integrator.
	/************************************************************************/
	STEP_RKF45,

	/************************************************************************/
	/* Embedded Runge-Kutta Cash-Karp (4, 5) method.
	/************************************************************************/
	STEP_CK,

	/************************************************************************/
	/* Embedded Runge-Kutta Prince-Dormand (8,9) method
	/************************************************************************/
	STEP_RK8PD,

	/************************************************************************/
	/* Implicit 2nd order Runge-Kutta at Gaussian points.
	/************************************************************************/
	STEP_RK2IMP,

	/************************************************************************/
	/* Implicit 4th order Runge-Kutta at Gaussian points.
	/************************************************************************/
	STEP_RK4IMP,

	/************************************************************************/
	/* Implicit Bulirsch-Stoer method of Bader and Deuflhard. This algorithm requires the Jacobian.
	/************************************************************************/
//	STEP_BSIMP,

	/************************************************************************/
	/* M=1 implicit Gear method.
	/************************************************************************/
	STEP_GEAR1,

	/************************************************************************/
	/* M=2 implicit Gear method.
	/************************************************************************/
//	STEP_GEAR2
}Enum_Step_Type;

typedef enum
{
	/************************************************************************/
	/* standard control
	/************************************************************************/
	STANDARD = 0,
	/************************************************************************/
	/* scaled control
	/************************************************************************/
	SCALED

}Enum_Control_Type;

#endif // NWPU_ODEIV_STEP_TYPE_H__