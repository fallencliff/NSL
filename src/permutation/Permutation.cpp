/********************************************************************
	filename: 	Permutation.cpp
	author:		huzhijian
	created:	5:5:2010   16:45
	brief:	≈≈¡–¿‡
*********************************************************************/

#include <gsl_errno.h>

#include <Permutation.h>

namespace gslcpp
{

	gsl_permutation* CPermutation::Alloc(const size_t n)
	{
		gsl_permutation * p;
		
		if (n == 0)
		{
			GSL_ERROR_VAL ("permutation length n must be positive integer",
				GSL_EDOM, 0);
		}
		
		p = (gsl_permutation *) new gsl_permutation;
		
		if (p == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for permutation struct",
				GSL_ENOMEM, 0);
		}
		
		p->data = (size_t *) new size_t[n];
		
		if (p->data == 0)
		{
			delete p;         /* exception in constructor, avoid memory leak */
			
			GSL_ERROR_VAL ("failed to allocate space for permutation data",
				GSL_ENOMEM, 0);
		}
		
		p->size = n;
		
		return p;
	}

	gsl_permutation * CPermutation::Calloc (const size_t n)
	{
		size_t i;
		
		gsl_permutation * p =  Alloc (n);
		
		if (p == 0)
			return 0;
		
		/* initialize permutation to identity */
		
		for (i = 0; i < n; i++)
		{
			p->data[i] = i;
		}
		
		return p;
	}

	CPermutation::CPermutation(size_t n): gsldata(NULL)
	{
		gsldata =Calloc(n);
		
	}

	CPermutation::CPermutation() : gsldata(NULL)
	{
		;
	}

	void CPermutation::Free()
	{
		if (gsldata)
		{
			if (gsldata->data)
			{
				delete [] gsldata->data;
				gsldata->data = NULL;
			}
			delete gsldata;
			gsldata = NULL;
		}
	}

	size_t CPermutation::Get(size_t i)
	{			
		if (i >= Size())         /* size_t is unsigned, can't be negative */
		{
			GSL_ERROR_VAL ("index out of range", GSL_EINVAL, 0);
		}			
		return gsldata->data[i];
	}

	void CPermutation::Resize(size_t n)
	{
		Free();
		gsldata= Calloc(n);
	}

	int CPermutation::Swap (const size_t i, const size_t j)
	{
		const size_t size = gsldata->size ;
		
		if (i >= size)
		{
			GSL_ERROR("first index is out of range", GSL_EINVAL);
		}
		
		if (j >= size)
		{
			GSL_ERROR("second index is out of range", GSL_EINVAL);
		}
		
		if (i != j)
		{
			size_t tmp = gsldata->data[i];
			gsldata->data[i] = gsldata->data[j];
			gsldata->data[j] = tmp;
		}
		
		return GSL_SUCCESS;
	}
	
	

	int CPermutation::Permute(double * data, const size_t stride, size_t n)
	{
		size_t i, k, pk;

		size_t* p = (gsldata->data);
		
		for (i = 0; i < n; i++)
		{
			k = p[i];
			
			while (k > i) 
				k = p[k];
			
			if (k < i)
				continue ;
			
			/* Now have k == i, i.e the least in its cycle */
			
			pk = p[k];
			
			if (pk == i)
				continue ;
			
			/* shuffle the elements of the cycle */
			
			{				
				double t = data[i*stride];
				
				while (pk != i)
				{
					
					double r1 = data[pk*stride];
					data[k*stride] = r1;
					
					k = pk;
					pk = p[k];
				};
								
				data[k*stride] = t;
			}
		}
		
		return GSL_SUCCESS;
	}

	void CPermutation::Identity()
	{
		const size_t n = Size();
		size_t* data = gsldata->data;
		/* initialize permutation to identity */
		
		for (size_t i = 0; i < n; i++)
		{
			data[i] = i;
		}
	}

	bool CPermutation::IsValid()
	{
		const size_t size = Size();
		const size_t* data = GetDataPtr();
		
		size_t i, j ;
		
		for (i = 0; i < size; i++) 
		{
			if (data[i] >= size)
			{
				//GSL_ERROR("permutation index outside range", GSL_FAILURE) ;
				return false;
			}
			
			for (j = 0; j < i; j++)
			{
				if (data[i] == data[j])
				{
					//GSL_ERROR("duplicate permutation index", GSL_FAILURE) ;
					return false;
				}
			}
		}
		
		return true;
	}

	void CPermutation::Reverse()
	{
		const size_t size = Size() ;
		size_t* data = GetDataPtr();

		for (size_t i = 0; i < (size / 2); i++) 
		{
			size_t j = size - i - 1;
			
			size_t tmp = data[i] ;
			data[i] = data[j] ;
			data[j] = tmp ;
		}
	}

	CPermutation& CPermutation::operator=(const CPermutation& other)
	{
		size_t size = other.Size();
		size_t* dest_data = GetDataPtr();
		const size_t* other_data = other.GetDataPtr();
		Resize(size);
		
		for (size_t j = 0; j < size; j++)
		{
			dest_data[j] = other_data[j];
		}

		return *this;
	}

	CPermutation CPermutation::Inverse()
	{
		const size_t size = Size();
				
		CPermutation inv(size);

		size_t* inv_data = inv.GetDataPtr();
		size_t* data = GetDataPtr();
		
		for (size_t i = 0; i < size; i++) 
		{
			inv_data[data[i]] = i ;
		}
		
		return inv;
	}

	int CPermutation::Next()
	{
	/* Replaces p with the next permutation (in the standard lexicographical
	* ordering).  Returns GSL_FAILURE if there is no next permutation.
		*/
		const size_t size = Size();
		size_t* data = GetDataPtr();
		size_t i, j, k;
		
		if (size < 2)
		{
			return GSL_FAILURE;
		}
		
		i = size - 2;
		
		while ((data[i] > data[i+1]) && (i != 0))
		{
			i--;
		}
		
		if ((i == 0) && (data[0] > data[1]))
		{
			return GSL_FAILURE;
		}
		
		k = i + 1;
		
		for (j = i + 2; j < size; j++ )
		{
			if ((data[j] > data[i]) && (data[j] < data[k]))
			{
				k = j;
			}
		}
		
		/* swap i and k */
		
		{
			size_t tmp = data[i];
			data[i] = data[k];
			data[k] = tmp;
		}
		
		for (j = i + 1; j <= ((size + i) / 2); j++)
		{
			size_t tmp = data[j];
			data[j] = data[size + i - j];
			data[size + i - j] = tmp;
		}
		
		return GSL_SUCCESS;
	}

	int CPermutation::Prev()
	{
		const size_t size = Size();
		size_t* data = GetDataPtr();
		size_t i, j, k;
		
		if (size < 2)
		{
			return GSL_FAILURE;
		}
		
		i = size - 2;
		
		while ((data[i] < data[i+1]) && (i != 0))
		{
			i--;
		}
		
		if ((i == 0) && (data[0] < data[1]))
		{
			return GSL_FAILURE;
		}
		
		k = i + 1;
		
		for (j = i + 2; j < size; j++ )
		{
			if ((data[j] < data[i]) && (data[j] > data[k]))
			{
				k = j;
			}
		}
		
		/* swap i and k */
		
		{
			size_t tmp = data[i];
			data[i] = data[k];
			data[k] = tmp;
		}
		
		for (j = i + 1; j <= ((size + i) / 2); j++)
		{
			size_t tmp = data[j];
			data[j] = data[size + i - j];
			data[size + i - j] = tmp;
		}
		
		return GSL_SUCCESS;
	}

	CPermutation CPermutation::Mul(const CPermutation& other)
	{
		size_t i;
		const size_t size_pa = Size();
		const size_t size_pb = other.Size();
		
		if (size_pa != size_pb)
		{
			GSL_ERROR("size of pb does not match size of pa", GSL_EINVAL);
		}
		
		CPermutation res(size_pa);
		size_t* res_data = res.GetDataPtr();
		const size_t* pa_data = GetDataPtr();
		const  size_t* pb_data = other.GetDataPtr();

		for (i = 0; i < size_pa; i++)
		{
			res_data[i] = pb_data[pa_data[i]];
		}
		
		return res;
	}

	CPermutation CPermutation::LinearToCanonical() const
	{
		const size_t n = Size();
		size_t i, k, s;
		size_t t = n;
		
		CPermutation q(n);

		const size_t *const pp = GetDataPtr();
		size_t *const qq = q.GetDataPtr();
		
			
		for (i = 0; i < n; i++)
		{
			
			k = pp[i];
			s = 1;
			
			while (k > i)
			{
				k = pp[k];
				s++;
			}
			
			if (k < i)
				continue;
			
			/* Now have k == i, i.e the least in its cycle, and s == cycle length */
			
			t -= s;
			
			qq[t] = i;
			
			k = pp[i];
			s = 1;
			
			while (k > i)
			{
				qq[t + s] = k;
				k = pp[k];
				s++;
			}
			
			if (t == 0)
				break;
		}
		
		return q;
	}

	CPermutation CPermutation::CanonicalToLinear() const
	{
		size_t i, k, kk, first;
		const size_t n = Size();
		
		CPermutation p(n);

		size_t *const pp = p.GetDataPtr();
		const size_t *const qq = GetDataPtr();
		
		for (i = 0; i < n; i++)
		{
			pp[i] = i;
		}
		
		k = qq[0];
		first = pp[k];
		
		for (i = 1; i < n; i++)
		{
			kk = qq[i];
			
			if (kk > first)
			{
				pp[k] = pp[kk];
				k = kk;
			}
			else
			{
				pp[k] = first;
				k = kk;
				first = pp[kk];
			}
		}
		
		pp[k] = first;
		
		return p;
	}

	size_t CPermutation::Inversions()
	{
		size_t count = 0;
		size_t i, j;
		const size_t size = Size();
		size_t* data = GetDataPtr();

		for (i = 0; i < size - 1; i++)
		{
			for (j = i + 1; j < size; j++)
			{
				if (data[i] > data[j])
				{
					count++;
				}
			}
		}
		
		return count;
	}

	size_t CPermutation::LinearCycles()
	{
		size_t i, k;
		size_t count = 0;
		const size_t size = Size();
		size_t* data = GetDataPtr();

		for (i = 0; i < size; i++)
		{
			
			k = data[i];
			
			while (k > i)
			{
				k = data[k];
			}
			
			if (k < i)
				continue;
			
			count++;
		}
		
		return count;
	}

	size_t CPermutation::CanonicalCycles()
	{
		size_t i;
		size_t count = 1;
		const size_t size = Size();
		size_t* data = GetDataPtr();

		size_t min = data[0];
		
		for (i = 0; i < size; i++)
		{
			if (data[i] < min)
			{
				min = data[i];
				count++;
			}
		}
		
		return count;
	}
}