/********************************************************************
	filename: 	RandGenerator.cpp
	author:		huzhijian
	created:	5:5:2010   16:55
	brief:	random number generator
*********************************************************************/
#include <iostream>
#include <algorithm>
#include <gsl_errno.h>
#include <RandGenerator.h>
#include <random/RandAlgorithm.h>
using std::cout;
using std::endl;
using std::transform;
namespace gslcpp
{

	const  CRandom::TypeArray CRandom::type_stack = CRandom::TypeInit(); //初始化生成器的类型堆栈
	
	CRandom::TypeArray CRandom::TypeInit()
	{
		CRandom::TypeArray res;
		res.push_back(gsl_rng_borosh13);
		res.push_back(gsl_rng_coveyou);
		res.push_back(gsl_rng_cmrg);
		res.push_back(gsl_rng_fishman18);
		res.push_back(gsl_rng_fishman20);
		res.push_back(gsl_rng_fishman2x);
		res.push_back(gsl_rng_gfsr4);
		res.push_back(gsl_rng_knuthran);
		res.push_back(gsl_rng_knuthran2);
		res.push_back(gsl_rng_lecuyer21);
		res.push_back(gsl_rng_minstd);
		res.push_back(gsl_rng_mrg);
		res.push_back(gsl_rng_mt19937);
		res.push_back(gsl_rng_mt19937_1999);
		res.push_back(gsl_rng_mt19937_1998);
		res.push_back(gsl_rng_r250);
		res.push_back(gsl_rng_ran0);
		res.push_back(gsl_rng_ran1);
		res.push_back(gsl_rng_ran2);
		res.push_back(gsl_rng_ran3);
		res.push_back(gsl_rng_rand);
		res.push_back(gsl_rng_rand48);
		res.push_back(gsl_rng_random128_bsd);
		res.push_back(gsl_rng_random128_glibc2);
		res.push_back(gsl_rng_random128_libc5);
		res.push_back(gsl_rng_random256_bsd);
		res.push_back(gsl_rng_random256_glibc2);
		res.push_back(gsl_rng_random256_libc5);
		res.push_back(gsl_rng_random32_bsd);
		res.push_back(gsl_rng_random32_glibc2);
		res.push_back(gsl_rng_random32_libc5);
		res.push_back(gsl_rng_random64_bsd);
		res.push_back(gsl_rng_random64_glibc2);
		res.push_back(gsl_rng_random64_libc5);
		res.push_back(gsl_rng_random8_bsd);
		res.push_back(gsl_rng_random8_glibc2);
		res.push_back(gsl_rng_random8_libc5);
		res.push_back(gsl_rng_random_bsd);
		res.push_back(gsl_rng_random_glibc2);
		res.push_back(gsl_rng_random_libc5);
		res.push_back(gsl_rng_randu);
		res.push_back(gsl_rng_ranf);
		res.push_back(gsl_rng_ranlux);
		res.push_back(gsl_rng_ranlux389);
		res.push_back(gsl_rng_ranlxd1);
		res.push_back(gsl_rng_ranlxd2);
		res.push_back(gsl_rng_ranlxs0);
		res.push_back(gsl_rng_ranlxs1);
		res.push_back(gsl_rng_ranlxs2);
		res.push_back(gsl_rng_ranmar);
		res.push_back(gsl_rng_slatec);
		res.push_back(gsl_rng_taus);
		res.push_back(gsl_rng_taus2);
		res.push_back(gsl_rng_taus113);
		res.push_back(gsl_rng_transputer);
		res.push_back(gsl_rng_tt800);
		res.push_back(gsl_rng_uni);
		res.push_back(gsl_rng_uni32);
		res.push_back(gsl_rng_vax);
		res.push_back(gsl_rng_waterman14);
		res.push_back(gsl_rng_zuf);

		return res;
	}
	CRandom::CRandom(RNG_ORDER type /*= RNG_TAUS*/)
	{
		if (type > type_stack.size() || type < RNG_ORDER(0))
		{
			gsl_error("type of generator not exists",
				__FILE__, __LINE__, GSL_ENOMEM);
		}

		/* FIXME: we're assuming that a char is 8 bits */
		state = (void *)new char[type_stack[type]->size];
	
		if (state == NULL)
		{		
			gsl_error("failed to allocate space for rng state",
				__FILE__, __LINE__, GSL_ENOMEM);
		}

		rng_type = type_stack[type];
		SetSeed(0);
	}

	CRandom::~CRandom()
	{
		Free();
		
	}
	void CRandom::PrintState()
	{
		char *p = (char *) (state);
		const size_t n = rng_type->size;
		
		for (size_t i = 0; i < n; i++)
		{
			/* FIXME: we're assuming that a char is 8 bits */
			printf ("%x\t", *(p + i));
		}
	}
	
	void CRandom::SetSeed(unsigned long int seed)
	{
		(rng_type->set)(state, seed);
	}
	
	unsigned long int CRandom::Rand(bool b_pos)
	{
		double x ;
		do
		{
			x = (rng_type->get)(state);
		}
		while (x == 0 && b_pos);
		
		return x ;
	}

	double CRandom::RandUniform(bool b_pos)
	{
		double x ;
		do
		{
			x = (rng_type->get_double)(state);
		}
		while (x == 0 && b_pos) ;
		
		return x ;
	}


	int CRandom::RandUniformInt(unsigned long int n, bool b_pos)
	{
		unsigned long int offset = Min();

		unsigned long int range = Max() - offset;
		
		
		if (n > range || n == 0) 
		{
			gsl_error("invalid n, either 0 or exceeds maximum value of generator",
				__FILE__, __LINE__, 0) ;
		}
		
		unsigned long int scale = range / n;
		unsigned long int k;

		do
		{
			do 
			{
				k = ( (rng_type->get)(state) - offset ) / scale;
			}
			while (k == 0 && b_pos) ;
			
		}
		while (k >= n);
		
		return k;
	}

	inline const char* CRandom::Name()
	{
		return rng_type->name;
	}
	
	inline unsigned long int CRandom::Max()
	{
		return rng_type->max;
	}
	
	inline unsigned long int CRandom::Min()
	{
		return rng_type->min;
	}

	inline size_t CRandom::Size()
	{
		return rng_type->size;
	}

	inline void* CRandom::State()
	{
		return state;
	}

	CRandom& CRandom::operator=(CRandom& other)
	{
		//string name = Name();
		//string other_name = other.Name();
		//if (name.compare(other_name) != 0)
		if (GetType() != other.GetType())
		{
			gsl_error ("generators must be of the same type",
                        __FILE__, __LINE__, GSL_EINVAL);
		}

		memcpy(this->state, other.State(), other.Size());

		return *this;
	}
	
	int CRandom::Save(char* filename)
	{
		FILE* fp = fopen(filename, "w");

		if (NULL == fp)
		{
			GSL_ERROR ("file open failed", GSL_EFAILED);
		}

		size_t n = rng_type->size ;
		
		char * m_state = (char *)state;
		
		size_t items = fwrite (m_state, 1, n, fp);
		
		if (items != n)
		{
			GSL_ERROR ("fwrite failed", GSL_EFAILED);
		}
		
		fclose(fp);
		return GSL_SUCCESS;
	}

	int CRandom::Load(string filename)
	{
		return GSL_SUCCESS;
	}
	
	const gsl_rng_type* CRandom::GetType()
	{
		return rng_type;	
	}
	
	bool CRandom::ChangeType( RNG_ORDER type )
	{
		Free();
		if (type > type_stack.size() || type < RNG_ORDER(0))
		{
			gsl_error("type of generator not exists",
				__FILE__, __LINE__, GSL_ENOMEM);
		}
		
		/* FIXME: we're assuming that a char is 8 bits */
		state = (void *)new char[type_stack[type]->size];
		
		if (state == NULL)
		{		
			gsl_error("failed to allocate space for rng state",
				__FILE__, __LINE__, GSL_ENOMEM);
		}
		
		rng_type = type_stack[type];

		return true;
	}
	
	void CRandom::Free()
	{
		if (state)
		{
			delete [] (char*)state;
			state = NULL;
		}	
	}
	
	void CRandom::printAllNames()
	{
		cout << "共有" <<type_stack.size() <<"种随机数算法" <<endl;
		size_t width_ = 25;
		cout.width(width_);
		cout << "名称";
		cout.width(width_);
		cout<<"枚举量" <<endl;
		for (size_t i=0; i<type_stack.size(); i++)
		{
			cout.width(width_);
			string name_ = type_stack[i]->name;
			cout << name_.c_str();
			int pos;
			if((pos = name_.find('-') )!= string::npos)
			{
				name_[pos] = '_';
			}
			transform(name_.begin(), name_.end(), name_.begin(), toupper);
			cout.width(width_);
			cout << "RNG_" << name_.c_str()<< endl;
		}
				
	}
}

