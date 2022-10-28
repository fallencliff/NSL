/********************************************************************
	filename: 	Ntuple.h
	author:		hu zhijian
	created:	5:5:2010   16:28
	brief:		Ntuple
*********************************************************************/
#ifndef NSL_NTUPLE_H__
#define NSL_NTUPLE_H__




#include <NSL.H>
#include <stdio.h>
#include <histogram.h>
#include <CommonStruct.h>

namespace gslcpp
{
	class NSL_EXPORT CNtuple
	{
	public:
		
	private:
	
		gsl_ntuple* gsldata;

		gsl_ntuple* Create(const char* filename, const char* mode, void* data, size_t size);
	public:
	
		CNtuple(const char* filename, const char* mode, void* data, size_t size);

		virtual ~CNtuple();

		int Free();

		int Write()const;
		
		int Read();

		int Project(CHistogram& h, value_fn* value_func, select_fn* select_func);

	};



}

#endif // NSL_NTUPLE_H__