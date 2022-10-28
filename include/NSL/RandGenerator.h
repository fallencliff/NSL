/********************************************************************
	filename: 	RandGenerator.h
	author:		hu zhijian
	created:	5:5:2010   16:55
	brief:	random number generator
*********************************************************************/

#ifndef NSL_RANDGENERATOR_H__
#define NSL_RANDGENERATOR_H__


#pragma warning( disable : 4251) //type_stack�о���


#include <NSL.h>
#include <vector>
#include <rgn_type.h>

using std::vector;
using std::string;


namespace gslcpp
{
	class  NSL_EXPORT CRandom
	{
	public:	

		typedef vector<const gsl_rng_type*> TypeArray;

		static const TypeArray type_stack;

		static TypeArray TypeInit();	

	private:
		void* state;
		const gsl_rng_type* rng_type;
		
		void Free();

	public:
		CRandom(RNG_ORDER type = RNG_TAUS);

		virtual ~CRandom();

		/*set seed*/
		void SetSeed(unsigned long int seed);
		
		/*����[min,max]�����һ������, b_posΪtrueʱ,���ط���ֵ*/
		unsigned long int Rand(bool b_pos = false);
	
		/*����[0,1)�����еĸ�����, b_posΪtrueʱ,���ط���ֵ*/
		double RandUniform(bool b_pos = false);
		
		/*����[0,n-1]�����е�����, b_posΪtrueʱ,���ط���ֵ*/
		int RandUniformInt(unsigned long int n, bool b_pos = false);
		
		const gsl_rng_type* GetType();

		const char* Name();	
		void printAllNames();
		unsigned long int Max();
		unsigned long int Min();
		size_t Size();

		void* State();
		void PrintState();	

		int Save(char* filename);
		int Load(string filename);

		CRandom& operator=(CRandom& other);

		bool ChangeType(RNG_ORDER type);

	};
}

#endif // NSL_RANDGENERATOR_H__
