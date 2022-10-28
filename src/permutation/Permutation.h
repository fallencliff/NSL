#ifndef NWPU_PERMUTATION_H
#define NWPU_PERMUTATION_H
/** 
 * @file Permutation.h.h
 * @brief 排列类  //缺少读写文件功能
 * @author 
 * @date 2010-03-31
 */

#include <SL_dll.h>

#include <Vector.h>
#include <Matrix.h>

namespace gslcpp
{
	
	
	class SL_DLL_API CPermutation
	{
	public:
		typedef struct
		{
			size_t size;
			size_t *data;
		}gsl_permutation;

	private:
		gsl_permutation* gsldata;

		void Free()
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
		
	public:
		CPermutation(size_t n)
		{
			gsldata =Calloc(n);

		}
		CPermutation() : gsldata(NULL)
		{
			;
		}
		~CPermutation()
		{
			Free();
			
		}
		
		void Resize(size_t n)
		{
			Free();
			gsldata= Calloc(n);
		}
		
		gsl_permutation* GslObj() const
		{
			return gsldata;
		}

		size_t Size() const
		{
			return gsldata->size;
		}
		
		gsl_permutation* Alloc(const size_t n);

		gsl_permutation * Calloc (const size_t n);

		/************************************************************************/
		/* This function exchanges the i-th and j-th elements of the permutation p                                                                     */
		/************************************************************************/
		int Swap (const size_t i, const size_t j);
		
		int Permute_Vector(CVector& v);
		
		int Permute(double * data, const size_t stride, size_t n);

		void Identity();

		/************************************************************************/
		/* This function returns the value of the i-th element of the permutation p. 
		If i lies outside the allowed range of 0 to n-1 then the error handler is invoked and 0 is returned.                                                                     */
		/************************************************************************/
		size_t Get(size_t i)
		{			
			if (i >= Size())         /* size_t is unsigned, can't be negative */
			{
				GSL_ERROR_VAL ("index out of range", GSL_EINVAL, 0);
			}			
			return gsldata->data[i];
		}

		size_t* GetDataPtr()
		{
			return gsldata->data;
		}

		const size_t* GetDataPtr() const
		{
			return gsldata->data;
		}

		/************************************************************************/
		/* This function checks that the permutation p is valid. 
		The n elements should contain each of the numbers 0 to n-1 once and only once.                                                                     */
		/************************************************************************/
		bool IsValid();

		/************************************************************************/
		/* This function reverses the elements of the permutation p.                                                                      */
		/************************************************************************/
		void Reverse();

		CPermutation Inverse();

		CPermutation& operator=(const CPermutation& other);

		/************************************************************************/
		/* This function advances the permutation p to the next permutation in lexicographic order and returns GSL_SUCCESS.
		If no further permutations are available it returns GSL_FAILURE and leaves p unmodified. 
		Starting with the identity permutation and repeatedly applying this function will iterate through all possible permutations of a given order.                                                                     */
		/************************************************************************/
		int Next();

		/************************************************************************/
		/* This function steps backwards from the permutation p to the previous permutation in lexicographic order, returning GSL_SUCCESS. 
		If no previous permutation is available it returns GSL_FAILURE and leaves p unmodified.                                                                     */
		/************************************************************************/
		int Prev();

		/************************************************************************/
		/* This function combines the two permutations pa and pb into a single permutation p, 
		where p = pa . pb. 
		The permutation p is equivalent to applying pb first and then pa.                                                                     */
		/************************************************************************/
		CPermutation Mul(const CPermutation& other);

		CPermutation LinearToCanonical() const;
		CPermutation CanonicalToLinear() const;
		
		/************************************************************************/
		/* This function counts the number of inversions in the permutation p. An inversion
		is any pair of elements that are not in order. For example, the permutation 2031
		has three inversions, corresponding to the pairs (2,0) (2,1) and (3,1). The identity
		permutation has no inversions.                                                                     */
		/************************************************************************/
		size_t Inversions();

		/************************************************************************/
		/* This function counts the number of cycles in the permutation p, given in linear form                                                                     */
		/************************************************************************/
		size_t LinearCycles();

		/************************************************************************/
		/* This function counts the number of cycles in the permutation q, given in canonical form                                                                     */
		/************************************************************************/
		size_t CanonicalCycles();
	};

	
}


#endif