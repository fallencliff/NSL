/********************************************************************
	filename: 	MonteCarlo.h
	author:		hu zhijian
	created:	6:5:2010   9:29
	brief:	Monte Carlo Integration
*********************************************************************/


#ifndef NSL_MONTECARLO_H__
#define NSL_MONTECARLO_H__

#include <NSL.h>
#include <stdlib.h>
#include <state.h>
#include <RandGenerator.h>

namespace gslcpp
{
	typedef enum {PLAIN = 1, MISER, VEGAS}ALG;

	class NSL_EXPORT CMonteCarlo
	{
	//public:
	
	private:
		MonteCarlo_function mc_func; //�����ֺ����ĵ�ַ
		void* mc_state;	//״̬��ַ
		ALG mc_alg;	//�㷨���� PLAIN or MISER or VEGAS
		CRandom rng;
		size_t dim;
		void* InitPlain(const size_t dim);
		void* InitMiser(const size_t dim);
		void* InitVegas(const size_t dim);
		
		double PlainIntegrate(double* xl, double* xu, const size_t dim,  size_t calls, double& abserr);
		double MiserIntegrate(double* xl, double* xu, const size_t dim,  size_t calls, double& abserr);
		double VegasIntegrate(double* xl, double* xu, const size_t dim,  size_t calls, double& abserr);

		double EstimateCorrmc
			(
			const double xl[], const double xu[],
			size_t dim, size_t calls,
			double &abserr,
			double sigma_l[], 
			double sigma_r[]
			);

		void PrintLimit(double xl[], double xu[], unsigned long dim);
		void PrintHead(unsigned long num_dim, unsigned long calls, unsigned int it_num, unsigned int bins, unsigned int boxes);
		void PrintRes
			(
			unsigned int itr, 
			double res, double err, 
			double cum_res, double cum_err,
			double chi_sq
			);
		
		void PrintDist(unsigned long dim);
		void PrintGrid(unsigned long dim);
		void RefineGrid();
		void ResizeGrid(unsigned int bins);
		void RandomPoint
			(
			double x[], int bin[], double *bin_vol,
			const int box[], const double xl[], const double xu[]
            );

		void AccumulateDistribution(int bin[], double y);
		void ResetGridValues();
		void InitGrid(const double xl[], const double xu[], size_t dim);
		int ChangeBoxCoord(int box[]);
		void InitBoxCoord(int box[]);
		/************************************************************************/
		/* ����Free�ͷ�state�ṹ����ָ��ָ����ڴ棬�����ͷ�stateָ����ڴ棬�������ͷ�
		/************************************************************************/
		void Free(void *state);
	public:
		CMonteCarlo(const size_t dimension, MonteCarlo_function F, RNG_ORDER rng_type = RNG_TAUS,ALG alg = PLAIN);
		virtual ~CMonteCarlo();

		void ResetFunc(MonteCarlo_function F);

		void ResetRngType(RNG_ORDER rng_type);

		void ResetAlg(ALG alg, const size_t dim);

		void* GetState() const;

		double Integrate(double* xl, double* xu,  const size_t calls, double& abserr);
	};
}


#endif // NSL_MONTECARLO_H__