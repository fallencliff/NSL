/********************************************************************
	filename: 	Combination.cpp
	author:		huzhijian
	created:	5:5:2010   16:50
	brief:	组合类
*********************************************************************/

#include <Combination.h>
#include <gsl_errno.h>
#include <stdio.h> //删除时和print函数一起删掉
namespace gslcpp
{
	CCombination::CCombination(const size_t n, const size_t k, bool flag)
	{
		gsldata = flag ? Alloc(n, k) : Calloc(n, k);
	}

	CCombination::~CCombination()
	{
		Free();
	}

	void CCombination::Free()
	{
		if (gsldata)
		{
			if (gsldata->k > 0 )
			{
				delete [] gsldata->data;
				gsldata->data = NULL;
			}
			delete gsldata;
			gsldata = NULL;
		}
	}

	gsl_combination* CCombination::Alloc(const size_t n, const size_t k)
	{
		gsl_combination * c;
		
		if (n == 0)
		{
			GSL_ERROR_VAL ("combination parameter n must be positive integer",
				GSL_EDOM, 0);
		}
		if (k > n)
		{
			GSL_ERROR_VAL ("combination length k must be an integer less than or equal to n",
				GSL_EDOM, 0);
		}
		c = (gsl_combination *) new gsl_combination;
		
		if (c == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for combination struct",
				GSL_ENOMEM, 0);
		}
		
		if (k > 0)
		{
			c->data = (size_t *) new size_t[k];
			
			if (c->data == 0)
			{
				delete c;             /* exception in constructor, avoid memory leak */
				
				GSL_ERROR_VAL ("failed to allocate space for combination data",
					GSL_ENOMEM, 0);
			}
		}
		else
		{
			c->data = 0;
		}
		
		c->n = n;
		c->k = k;
		
		return c;
	}

	gsl_combination* CCombination::Calloc(const size_t n, const size_t k)
	{
		
		gsl_combination * c =  Alloc (n, k);
		
		if (c == 0)
			return 0;
		
		/* initialize combination to identity */
		
		for (size_t i = 0; i < k; i++)
		{
			c->data[i] = i;
		}
		
		return c;
	}
	
	void CCombination::InitFirst()
	{
		const size_t k = gsldata->k ;
		size_t* data = gsldata->data;
		/* initialize combination to identity */
		
		for (size_t i= 0; i < k; i++)
		{
			data[i] = i;
		}
		
	}
	void CCombination::InitLast()
	{
		const size_t k = gsldata->k ;
		size_t n = gsldata->n;
		size_t* data = gsldata->data;
		/* initialize combination to identity */
		
		for (size_t i = 0; i < k; i++)
		{
			data[i] = n - k + i;
		}
	}

	void CCombination::Resize(size_t n, size_t k, bool flag)
	{
		Free();
	
		gsldata = flag ? Alloc(n, k) : Calloc(n, k);
	}

	CCombination& CCombination::operator=(const CCombination& other)
	{
		
		const size_t src_n = other.gsldata->n;
		const size_t src_k = other.gsldata->k;
		const size_t dest_n = gsldata->n;
		const size_t dest_k = gsldata->k;
		
		if (src_n != dest_n || src_k != dest_k)
		{
			Resize(src_n, src_k);
		}
		
		const size_t* src_data = other.gsldata->data;
		size_t* dest_data = gsldata->data;

		for (size_t j = 0; j < src_k; j++)
		{
			dest_data[j] = src_data[j];
		}
		
		return *this;
	}

	inline size_t CCombination::operator[](const size_t index) const
	{
		
		if (index >= gsldata->k)            /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL ("index out of range", GSL_EINVAL, 0);
		}
		
		
		return gsldata->data[index];
	}

	inline bool CCombination::IsValid()
	{
		const size_t n = gsldata->n ;
		const size_t k = gsldata->k ;
		const size_t* data = gsldata->data;
		
		if( k > n )
		{
			return false;
		}
		for (size_t i = 0; i < k; i++) 
		{
			const size_t ci = data[i];
			
			if (ci >= n)
			{
				return false;
			}
			
			for (size_t j = 0; j < i; j++)
			{
				if (data[j] == ci)
				{
					return false;
				}
				if (data[j] > ci)
				{
					return false;
				}
			}
		}
		
		return true;
	}

	bool CCombination::Next()
	{
	/* Replaces c with the next combination (in the standard lexicographical
	* ordering).  Returns GSL_FAILURE if there is no next combination.
		*/
		const size_t n = gsldata->n;
		const size_t k = gsldata->k;
		size_t *data = gsldata->data;
		size_t i;
		
		if(k == 0)
		{
			return false;
		}
		i = k - 1;
		
		while(i > 0 && data[i] == n - k + i)
		{
			i--;
		}
		if(i == 0 && data[i] == n - k)
		{
			return false;
		}
		data[i]++;
		for(; i < k - 1; i++)
		{
			data[i + 1] = data[i] + 1;
		}
		return true;

	}
	bool CCombination::Prev()
	{
	/* Replaces c with the previous combination (in the standard
	* lexicographical ordering).  Returns GSL_FAILURE if there is no
	* previous combination.
		*/
		const size_t n = gsldata->n;
		const size_t k = gsldata->k;
		size_t *data = gsldata->data;
		size_t i;
		
		if(k == 0)
		{
			return false;
		}
		i = k - 1;
		
		while(i > 0 && data[i] == data[i-1] + 1)
		{
			i--;
		}
		if(i == 0 && data[i] == 0)
		{
			return false;
		}
		data[i++]--;
		for(; i < k; i++)
		{
			data[i] = n - k + i;
		}
		return true;

	}

	bool CCombination::Save(const char* filename, const char* format) const
	{
		int status;
		
		FILE* fp = fopen(filename, "w");
		
		size_t size = gsldata->k;
		
		const size_t* data = gsldata->data;
		
		for (size_t i=0; i<size; i++)
		{
			
			status = fprintf(fp, format, data[i]);
			if (status < 0)
			{
				gsl_error("fprintf failed", __FILE__, __LINE__, GSL_EFAILED);
			}
			fprintf(fp, "  ");
		}
		
		fprintf(fp, "\n");
		
		fclose(fp);
		
		return true;
	}

	bool CCombination::Load(const char* filename)
	{
		size_t k = gsldata->k ;
		
		size_t * data = gsldata->data ;
		
		FILE* fp = fopen(filename, "w");

		for (size_t i = 0; i < k; i++)
		{
			unsigned long j ;  
			
			/* FIXME: what if size_t != unsigned long ??? 
			
			  want read in size_t but have to read in unsigned long to avoid
			error from compiler */
			
			int status = fscanf (fp, "%u", &j);  
			
			if (status != 1)
			{
				gsl_error("fscanf failed", __FILE__, __LINE__, GSL_EFAILED);
			}
			
			data[i] = j;
		}
		
		fclose(fp);

		return true;
	}

	void CCombination::print()
	{
		size_t size = gsldata->k;
		
		const size_t* data = gsldata->data;
		
		for (size_t i=0; i<size; i++)
		{
			
			printf("%u ", data[i]);
		}
		printf("\n");
	}
	
}