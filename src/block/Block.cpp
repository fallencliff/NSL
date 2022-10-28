/********************************************************************
	filename: 	Block.cpp
	author:		hu zhijian
	created:	12:5:2010   20:09
	brief:	a class for vector and matrix
*********************************************************************/


#include <gsl_errno.h>
#include <Block.h>


namespace gslcpp
{
	CBlock::CBlock(const size_t n, bool clear): gsldata(NULL)
	{
		if (clear)
		{		
			gsldata = calloc(n);		
		} 
		else
		{
			gsldata = alloc(n);
		}
	}
	
	CBlock::CBlock(const CBlock& other): gsldata(NULL)
	{
		this->gsldata = other.gsldata;
	}

	CBlock::~CBlock()
	{
		if(gsldata)
		{
			delete [] gsldata->data;
			delete gsldata;
		}
		
	}

	size_t CBlock::Size() const
	{
		return this->gsldata->size;
	}

	gsl_block* CBlock::alloc(const size_t n)
	{
		gsl_block * b;
		
		if (n == 0)
		{
			GSL_ERROR_VAL ("block length n must be positive integer",
				GSL_EINVAL, 0);
		}
		
 		b = new(gsl_block);
		
		if (b == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for block struct",
				GSL_ENOMEM, 0);
		}
		
		b->data = (double *) new double[n];
		
		if (b->data == 0)
		{
			delete b;         /* exception in constructor, avoid memory leak */
 			
			GSL_ERROR_VAL ("failed to allocate space for block data",
				GSL_ENOMEM, 0);
		}
		
		b->size = n;

		return b;
	}

	gsl_block* CBlock::calloc(const size_t n)
	{
		size_t i;
		
		gsl_block * b = alloc(n);
		
		if (b == 0)
			return 0;
		
		/* initialize block to zero */
		
		for (i = 0; i < n; i++)
		{
			b->data[i] = 0;
		}

		return b;
	}

	int CBlock::Save(const char* filename, const char* format) const
	{

		FILE * stream =fopen(filename, "w");
			
		size_t i;
		
		for (i = 0; i < gsldata->size; i++)
		{

			int status;
							
			status = fprintf (stream, format, gsldata->data[i]);

			if (status < 0)
			{
				GSL_ERROR ("fprintf failed", GSL_EFAILED);
			}

			fprintf(stream, " ");
			
		}
		
		fprintf(stream, "\n");
		
		fclose(stream);

		return 0;
		
	}

	int CBlock::RawSave(const char* filename, 
		const size_t n,
		const size_t stride, 
		const char *format)
	{
		
		FILE * stream =fopen(filename, "w");
		
		size_t i;
		
		for (i = 0; i < n; i++)
		{
			
			int status;
			
			status = fprintf (stream, format, gsldata->data[i*stride]);
			
			if (status < 0)
			{
				GSL_ERROR ("fprintf failed", GSL_EFAILED);
			}
			
			fprintf(stream, " ");
			
		}
		
		fprintf(stream, "\n");
		
		fclose(stream);
		
		return 0;
	}
	int CBlock::Load(const char* filename)
	{	

		FILE* fp = fopen(filename, "r");

		size_t i;
		
		for (i = 0; i < gsldata->size; i++)
		{		
			int status = fscanf (fp, "%lf", &gsldata->data[i]) ;

			if (status != 1)
			{
				GSL_ERROR ("fscanf failed", GSL_EFAILED);
			}
			
		}

		fclose(fp);

		return GSL_SUCCESS;
	}

	double * CBlock::GetDataPtr() const
	{
		return gsldata->data;
	}
	
}