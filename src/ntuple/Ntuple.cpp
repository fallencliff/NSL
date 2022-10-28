/********************************************************************
	filename: 	Ntuple.cpp
	author:		hu zhijian
	created:	5:5:2010   16:40
	brief:	
*********************************************************************/
#include <Ntuple.h>
#include <gsl_errno.h>

namespace gslcpp
{
	CNtuple::CNtuple(const char* filename, const char* mode, void* data, size_t size)
		:gsldata(NULL)
	{
		gsldata = Create(filename, mode, data, size);
	}

	gsl_ntuple* CNtuple::Create(const char* filename, const char* mode, void* data, size_t size)
	{
		gsl_ntuple *ntuple = (gsl_ntuple *) new gsl_ntuple;
		
		if (ntuple == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for ntuple struct",
				GSL_ENOMEM, 0);
		}
		
		ntuple->ntuple_data = data;
		ntuple->size = size;
		
		ntuple->file = fopen (filename, mode);
		
		if (ntuple->file == 0)
		{
			delete ntuple;
			GSL_ERROR_VAL ("unable to create ntuple file", GSL_EFAILED, 0);
		}
		
		return ntuple;
	}

	int CNtuple::Write() const
	{
		size_t nwrite;
		
		nwrite = fwrite (gsldata->ntuple_data, gsldata->size,
			1, gsldata->file);
		
		if (nwrite != 1)
		{
			GSL_ERROR ("failed to write ntuple entry to file", GSL_EFAILED);
		}
		
		return GSL_SUCCESS;
	}

	int CNtuple::Read()
	{
		size_t nread;
		
		nread = fread (gsldata->ntuple_data, gsldata->size, 1, gsldata->file);
		
		if (nread == 0 && feof(gsldata->file))
		{
			return GSL_EOF;
		}
		
		if (nread != 1)
		{
			GSL_ERROR ("failed to read ntuple entry from file", GSL_EFAILED);
		}
		
		return GSL_SUCCESS;
	}
	
	int CNtuple::Free()
	{
		int status = fclose (gsldata->file);
		
		if (status)
		{
			GSL_ERROR ("failed to close ntuple file", GSL_EFAILED);
		}
		
		delete gsldata;
		
		return GSL_SUCCESS;	
	}
	
	CNtuple::~CNtuple()
	{
		if (gsldata)
		{
			Free();
		}
		
	}
	
	int CNtuple::Project( CHistogram& h, value_fn* value_func, select_fn* select_func)
	{
		size_t nread;
		#define EVAL(f,x) ((*((f)->function))(x,(f)->params))
		do
		{
			nread = fread (gsldata->ntuple_data, gsldata->size,
				1, gsldata->file);
			
			if (nread == 0 && feof(gsldata->file))
			{
				break ;
			}
			
			if (nread != 1) 
			{
				GSL_ERROR ("failed to read ntuple for projection", GSL_EFAILED);
			}
			
			if (EVAL(select_func, gsldata->ntuple_data))
			{
				h.Increment(EVAL(value_func, gsldata->ntuple_data)) ;
			}
		}
		while (1);
		#undef EVAL
		return GSL_SUCCESS;	
	}
}