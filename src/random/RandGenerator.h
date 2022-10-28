/**
 *@file RandGenerator.h
 *@brief
 *@author
 *@date 2010-03-24
 */
#ifndef NWPU_RANDGENERATOR_H
#define NWPU_RANDGENERATOR_H

#include <SL_dll.h>

#include <cstring>
#include <vector>
#include <rgn_type.h>

using std::vector;
using std::string;

#pragma warning( disable : 4251) //type_stack�о���
namespace gslcpp
{
	class  SL_DLL_API CRandom
	{
	public:
		typedef struct
		{
			const char *name;
			unsigned long int max;
			unsigned long int min;
			size_t size;
			void (*set) (void *state, unsigned long int seed);
			unsigned long int (*get) (void *state);
			double (*get_double) (void *state);

		}gsl_rng_type;

		typedef vector<const gsl_rng_type*> TypeArray;

		static const TypeArray type_stack;

		static CRandom::TypeArray TypeInit();	

	private:
		void* state;
		const gsl_rng_type* rng_type;
		
	public:
		CRandom(RNG_ORDER type = RNG_TAUS);

		~CRandom();

		/*set seed*/
		void Set(unsigned long int seed);
		
		/*����[min,max]�����һ������, b_posΪtrueʱ,���ط���ֵ*/
		unsigned long int Rand(bool b_pos = false);
	
		/*����[0,1)�����еĸ�����, b_posΪtrueʱ,���ط���ֵ*/
		double RandUniform(bool b_pos = false);
		
		/*����[0,n-1]�����е�����, b_posΪtrueʱ,���ط���ֵ*/
		int RandUniformInt(unsigned long int n, bool b_pos = false);
		

		string Name();	
		unsigned long int Max();
		unsigned long int Min();
		size_t Size();
		void* State();
		void PrintState();	

		int Save(char* filename);
		int Load(string filename);

		CRandom& operator=(CRandom& other);
	};
}

#endif
