/********************************************************************
	filename: 	ode_include.h
	author:		hu zhijian
	created:	7:5:2010   14:45
	brief:	this header file include all the  algorithm file for ode intergration
*********************************************************************/
#ifndef NWPU_ODE_INCLUDE_H__
#define NWPU_ODE_INCLUDE_H__

//for step
#include <ode-initval/ode_alg/rk2.cpp>
#include <ode-initval/ode_alg/rk4.cpp>
#include <ode-initval/ode_alg/rkf45.cpp>
#include <ode-initval/ode_alg/rkck.cpp>
#include <ode-initval/ode_alg/rk8pd.cpp>
#include <ode-initval/ode_alg/rk2imp.cpp>
#include <ode-initval/ode_alg/rk4imp.cpp>
//#include <ode-initval/ode_alg/bsimp.cpp>
#include <ode-initval/ode_alg/gear1.cpp>
//#include <ode-initval/ode_alg/gear2.cpp>


//for control
#include <ode-initval/ode_alg/cstd.cpp>
#include <ode-initval/ode_alg/cscal.cpp>

#endif // NWPU_ODE_INCLUDE_H__
