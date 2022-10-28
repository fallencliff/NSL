#ifndef NWPU_RGN_TYPE_H
#define NWPU_RGN_TYPE_H
namespace gslcpp
{
	
	
	enum RNG_ORDER_T
	{
			RNG_BOROSH13,
			RNG_COVEYOU,
			RNG_CMRG,
			RNG_FISHMAN18,
			RNG_FISHMAN20,
			RNG_FISHMAN2X,
			RNG_GFSR4,
			RNG_KNUTHRAN,
			RNG_KNUTHRAN2,
			RNG_LECUYER21,
			RNG_MINSTD,
			RNG_MRG,
			RNG_MT19937,
			RNG_MT19937_1999,
			RNG_MT19937_1998,
			RNG_R250,
			RNG_RAN0,
			RNG_RAN1,
			RNG_RAN2,
			RNG_RAN3,
			RNG_RAND,
			RNG_RAND48,
			RNG_RANDOM128_BSD,
			RNG_RANDOM128_GLIBC2,
			RNG_RANDOM128_LIBC5,
			RNG_RANDOM256_BSD,
			RNG_RANDOM256_GLIBC2,
			RNG_RANDOM256_LIBC5,
			RNG_RANDOM32_BSD,
			RNG_RANDOM32_GLIBC2,
			RNG_RANDOM32_LIBC5,
			RNG_RANDOM64_BSD,
			RNG_RANDOM64_GLIBC2,
			RNG_RANDOM64_LIBC5,
			RNG_RANDOM8_BSD,
			RNG_RANDOM8_GLIBC2,
			RNG_RANDOM8_LIBC5,
			RNG_RANDOM_BSD,
			RNG_RANDOM_GLIBC2,
			RNG_RANDOM_LIBC5,
			RNG_RANDU,
			RNG_RANF,
			RNG_RANLUX,
			RNG_RANLUX389,
			RNG_RANLXD1,
			RNG_RANLXD2,
			RNG_RANLXS0,
			RNG_RANLXS1,
			RNG_RANLXS2,
			RNG_RANMAR,
			RNG_SLATEC,
			RNG_TAUS,
			RNG_TAUS2,
			RNG_TAUS113,
			RNG_TRANSPUTER,
			RNG_TT800,
			RNG_UNI,
			RNG_UNI32,
			RNG_VAX,
			RNG_WATERMAN14,
			RNG_ZUF,
	};
	
	typedef enum RNG_ORDER_T RNG_ORDER;
	
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

}

#endif