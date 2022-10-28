/********************************************************************
	filename: 	MonteCarlo.cpp
	author:		hu zhijian
	created:	6:5:2010   9:26
	brief:	Monte Carlo Integration
*********************************************************************/
#include <math.h>
#include <MathConst.h>
#include <MonteCarlo.h>
#include <gsl_errno.h>

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

namespace gslcpp
{
#define MONTE_FN_EVAL(F,x) (*((F).f))(x, dim,(F).params)

/* A separable grid with coordinates and values */
#define COORD(s,i,j) ((s)->xi[(i)*(s)->dim + (j)])
#define NEW_COORD(s,i) ((s)->xin[(i)])
#define VALUE(s,i,j) ((s)->d[(i)*(s)->dim + (j)])


	CMonteCarlo::CMonteCarlo(const size_t dimension, MonteCarlo_function F, 
							 RNG_ORDER rng_type/* = RNG_TAUS*/,  ALG alg /* = PLAIN*/)
		:mc_func(F), mc_alg(alg), mc_state(NULL), rng(rng_type), dim(dimension)
	{
		ResetAlg(alg, dim);
	}
	
	CMonteCarlo::~CMonteCarlo()
	{
		if(mc_state)
		{
			Free(mc_state);
			delete mc_state;
			mc_state = NULL;
		}
		
	}
	
	void* CMonteCarlo::InitPlain( const size_t dim )
	{
		monte_plain_state *s = new monte_plain_state;
		
		if (s == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for plain state struct",
				GSL_ENOMEM, 0);
		}
		
		s->x = new double[dim];
		
		if (s->x == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for working vector",
				GSL_ENOMEM, 0);
		}
		
		s->dim = dim;
		
		return (void*)s;
	}
	
	void* CMonteCarlo::InitMiser( const size_t dim )
	{
		monte_miser_state *s =new monte_miser_state;
		
		if (s == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for miser state struct",
				GSL_ENOMEM, 0);
		}
		
		s->x = new double[dim];
		
		if (s->x == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
		}
		
		s->xmid = new double[dim];
		
		if (s->xmid == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for xmid", GSL_ENOMEM, 0);
		}
		
		s->sigma_l = new double[dim];
		
		if (s->sigma_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for sigma_l", GSL_ENOMEM, 0);
		}
		
		s->sigma_r = new double[dim];
		
		if (s->sigma_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for sigma_r", GSL_ENOMEM, 0);
		}
		
		s->fmax_l = new double[dim];
		
		if (s->fmax_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fmax_l", GSL_ENOMEM, 0);
		}
		
		s->fmax_r = new double[dim];
		
		if (s->fmax_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fmax_r", GSL_ENOMEM, 0);
		}
		
		s->fmin_l = new double[dim];
		
		if (s->fmin_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fmin_l", GSL_ENOMEM, 0);
		}
		
		s->fmin_r = new double[dim];
		
		if (s->fmin_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fmin_r", GSL_ENOMEM, 0);
		}
		
		s->fsum_l = new double[dim];
		
		if (s->fsum_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum_l", GSL_ENOMEM, 0);
		}
		
		s->fsum_r = new double[dim];
		
		if (s->fsum_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum_r", GSL_ENOMEM, 0);
		}
		
		s->fsum2_l = new double[dim];
		
		if (s->fsum2_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum2_l", GSL_ENOMEM, 0);
		}
		
		s->fsum2_r = new double[dim];
		
		if (s->fsum2_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum2_r", GSL_ENOMEM, 0);
		}
		
		
		s->hits_r = new size_t[dim];
		
		if (s->hits_r == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum2_r", GSL_ENOMEM, 0);
		}
		
		s->hits_l = new size_t[dim];
		
		if (s->hits_l == 0)
		{
			Free( (void*) s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for fsum2_r", GSL_ENOMEM, 0);
		}
		
		s->dim = dim;
		
		/* We use 8 points in each halfspace to estimate the variance. There are
		2*dim halfspaces. A variance estimate requires a minimum of 2 points. */
		s->min_calls = 16 * s->dim;
		s->min_calls_per_bisection = 32 * s->min_calls;
		s->estimate_frac = 0.1;
		s->alpha = 2.0;
		s->dither = 0.0;
		
		return (void*)s;	
	}
	
	void CMonteCarlo::Free(void *m_state)
	{
		//void* state = 
		switch(mc_alg)
		{
		case PLAIN:
			{
				monte_plain_state* state =  (monte_plain_state*)m_state;
				if (state)
				{
					if (state->x)
					{
						delete [] state->x;
					}
					
				}
				break;
			}
		case MISER:
			{
				monte_miser_state* state = (monte_miser_state*) m_state;
				if (state)
				{
					if (state->x)
					{
						delete [] state->x;
					}
					if (state->xmid)
					{
						delete [] state->xmid;
					}
					if (state->sigma_l)
					{
						delete [] state->sigma_l;
					}
					if (state->sigma_r)
					{
						delete [] state->sigma_r;
					}
					if (state->fmax_l)
					{
						delete [] state->fmax_l;
					}
					if (state->fmax_r)
					{
						delete [] state->fmax_r;
					}
					if (state->fmin_l)
					{
						delete [] state->fmin_l;
					}
					if (state->fmin_r)
					{
						delete [] state->fmin_r;
					}
					if (state->fsum_l)
					{
						delete [] state->fsum_l;
					}
					if (state->fsum_r)
					{
						delete [] state->fsum_r;
					}
					if (state->fsum2_l)
					{
						delete [] state->fsum2_l;
					}
					if (state->fsum2_r)
					{
						delete [] state->fsum2_r;
					}
					if (state->hits_r)
					{
						delete [] state->hits_r;
					}
					if (state->hits_l)
					{
						delete [] state->hits_l;
					}
				}
				break;
			}
		case VEGAS:
			{
				monte_vegas_state* state = (monte_vegas_state*) m_state;
				if (state)
				{

					if (state->delx)
					{
						delete [] state->delx;
					}
					if (state->d)
					{
						delete [] state->d;
					}
					if (state->xi)
					{
						delete [] state->xi;
					}
					if (state->xin)
					{
						delete [] state->xin;
					}
					if (state->weight)
					{
						delete [] state->weight;
					}
					if (state->box)
					{
						delete [] state->box;
					}
					if (state->bin)
					{
						delete [] state->bin;
					}
					if (state->x)
					{
						delete [] state->x;
					}
				}
				break;
			}
		}
	}
	
	void* CMonteCarlo::InitVegas( const size_t dim )
	{
		#define BINS_MAX 50             /* even integer, will be divided by two */
		typedef int coord;

		monte_vegas_state *s = new  monte_vegas_state;
			
		if (s == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for vegas state struct",
				GSL_ENOMEM, 0);
		}
		
		s->delx = new double[dim];
		
		if (s->delx == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for delx", GSL_ENOMEM, 0);
		}
		
		s->d = new double[BINS_MAX * dim];
		
		if (s->d == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for d", GSL_ENOMEM, 0);
		}
		
		s->xi = new double[(BINS_MAX+1) * dim];
		
		if (s->xi == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for xi", GSL_ENOMEM, 0);
		}
		
		s->xin = new double[(BINS_MAX+1) * dim];
		
		if (s->xin == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for xin", GSL_ENOMEM, 0);
		}
		
		s->weight = new double[(BINS_MAX) * dim];
		
		if (s->weight == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for xin", GSL_ENOMEM, 0);
		}
		
		s->box = new coord[dim];
		
		if (s->box == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for box", GSL_ENOMEM, 0);
		}
		
		s->bin = new coord[dim];
		
		if (s->bin == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for bin", GSL_ENOMEM, 0);
		}
		
		s->x = new double[dim];
		
		if (s->x == 0)
		{
			Free( (void*)s);
			delete s;
			s = NULL;
			GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
		}
		
		s->dim = dim;
		s->bins_max = BINS_MAX;

		/* Set some default values and whatever */
		s->stage = 0;
		s->alpha = 1.5;
		s->verbose = -1;
		s->iterations = 5;
		s->mode = VEGAS_MODE_IMPORTANCE;
		s->chisq = 0;
		s->bins = s->bins_max;
		s->ostream = stdout;
		
		return s;
		#undef BINS_MAX
	}
	
	inline void CMonteCarlo::ResetFunc(MonteCarlo_function F)
	{
		mc_func = F;	
	}
	
	inline void CMonteCarlo::ResetRngType( RNG_ORDER rng_type )
	{
		rng.ChangeType(rng_type);	
	}
	
	void CMonteCarlo::ResetAlg(ALG alg , const size_t dim)
	{
		if (mc_state)
		{
			Free(mc_state);
			delete mc_state;
			mc_state = NULL;
		}
		this->dim = dim;
		this->mc_alg = alg;
		switch(alg)
		{
		case PLAIN:
			{
				mc_state = InitPlain(dim);
				break;
			}
		case MISER:
			{
				mc_state = InitMiser(dim);
				break;
			}
		case VEGAS:
			{
				mc_state = InitVegas(dim);
				break;
			}
		default:
			gsl_error("no MonteCarlo algorithm like this exists", __FILE__, __LINE__, GSL_EINVAL);
		}

	}
	
	double CMonteCarlo::Integrate( double* xl, double* xu, const size_t calls, double &abserr )
	{
		double res = 0.0;

		switch(mc_alg)
		{
		case PLAIN:
			{
				res = PlainIntegrate(xl, xu, dim, calls, abserr);
				break;
			}
		case MISER:
			{
				res = MiserIntegrate(xl, xu, dim, calls, abserr);
				break;
			}
		case VEGAS:
			{
				res = VegasIntegrate(xl, xu, dim, calls, abserr);
				break;
			}
		}

		return res;
	}
	
	double CMonteCarlo::PlainIntegrate( double* xl, double* xu, const size_t dim, size_t calls, double& abserr )
	{
		double vol, m = 0, q = 0;
		monte_plain_state* state = (monte_plain_state*)mc_state;
		double *x = state->x;
		size_t n, i;
		
		if (dim != state->dim)
		{
			GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
		}
		
		for (i = 0; i < dim; i++)
		{
			if (xu[i] <= xl[i])
			{
				GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
			}
			
			if (xu[i] - xl[i] >DOUBLE_MAX)
			{
				GSL_ERROR ("Range of integration is too large, please rescale",
					GSL_EINVAL);
			}
		}
		
		/* Compute the volume of the region */
		
		vol = 1;
		
		for (i = 0; i < dim; i++)
		{
			vol *= xu[i] - xl[i];
		}
		
		for (n = 0; n < calls; n++)
		{
			/* Choose a random point in the integration region */
			
			for (i = 0; i < dim; i++)
			{
				x[i] = xl[i] + rng.RandUniform(true) * (xu[i] - xl[i]);
			}
			
			{
				double fval = MONTE_FN_EVAL (mc_func, x);
				
				/* recurrence for mean and variance */
				
				double d = fval - m;
				m += d / (n + 1.0);
				q += d * d * (n / (n + 1.0));
			}
		}
		
		if (calls < 2)
		{
			abserr = CMath::POSINF;
		}
		else
		{
			abserr = vol * sqrt (q / (calls * (calls - 1.0)));
		}
		return  vol * m;
	}
	
	double CMonteCarlo::MiserIntegrate( double* xl, double* xu, const size_t dim,  size_t calls, double& abserr )
	{
		size_t n, estimate_calls, calls_l, calls_r;
		monte_miser_state* state = (monte_miser_state*)mc_state;
		const size_t min_calls = state->min_calls;
		size_t i;
		size_t i_bisect;
		int found_best;
		
		double res_est = 0, err_est = 0;
		double res_r = 0, err_r = 0, res_l = 0, err_l = 0;
		double xbi_l, xbi_m, xbi_r, s;
		
		double vol;
		double weight_l, weight_r;
		
		double *x = state->x;
		double *xmid = state->xmid;
		double *sigma_l = state->sigma_l, *sigma_r = state->sigma_r;
		
		if (dim != state->dim)
		{
			GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
		}
		
		for (i = 0; i < dim; i++)
		{
			if (xu[i] <= xl[i])
			{
				GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
			}
			
			if (xu[i] - xl[i] > DOUBLE_MAX)
			{
				GSL_ERROR ("Range of integration is too large, please rescale",
					GSL_EINVAL);
			}
		}
		
		if (state->alpha < 0)
		{
			GSL_ERROR ("alpha must be non-negative", GSL_EINVAL);
		}
		
		/* Compute volume */
		
		vol = 1;
		
		for (i = 0; i < dim; i++)
		{
			vol *= xu[i] - xl[i];
		}
		
		if (calls < state->min_calls_per_bisection)
		{
			double m = 0.0, q = 0.0;
			
			if (calls < 2)
			{
				GSL_ERROR ("insufficient calls for subvolume", GSL_EFAILED);
			}
			
			for (n = 0; n < calls; n++)
			{
				/* Choose a random point in the integration region */
				
				for (i = 0; i < dim; i++)
				{
					x[i] = xl[i] + rng.RandUniform(true) * (xu[i] - xl[i]);
				}
				
				{
					double fval = MONTE_FN_EVAL (mc_func, x);
					
					/* recurrence for mean and variance */
					
					double d = fval - m;
					m += d / (n + 1.0);
					q += d * d * (n / (n + 1.0));
				}
			}
			
			abserr = vol * sqrt (q / (calls * (calls - 1.0)));
			return  vol * m;
		}
		
		estimate_calls = MAX(min_calls, calls * (state->estimate_frac));
		
		if (estimate_calls < 4 * dim)
		{
			GSL_ERROR ("insufficient calls to sample all halfspaces", GSL_ESANITY);
		}
		
		/* Flip coins to bisect the integration region with some fuzz */
		
		for (i = 0; i < dim; i++)
		{
			s = (rng.RandUniform(true) - 0.5) >= 0.0 ? state->dither : -state->dither;
			state->xmid[i] = (0.5 + s) * xl[i] + (0.5 - s) * xu[i];
		}
		
		/* The idea is to chose the direction to bisect based on which will
		give the smallest total variance.  We could (and may do so later)
		use MC to compute these variances.  But the NR guys simply estimate
		the variances by finding the min and max function values 
		for each half-region for each bisection. */
		
 		res_est = EstimateCorrmc( xl, xu, dim, estimate_calls, err_est, sigma_l, sigma_r);
		
		/* We have now used up some calls for the estimation */
		
		calls -= estimate_calls;
		
		/* Now find direction with the smallest total "variance" */
		
		{
			double best_var = DOUBLE_MAX;
			double beta = 2.0 / (1.0 + state->alpha);
			found_best = 0;
			i_bisect = 0;
			weight_l = weight_r = 1.0;
			
			for (i = 0; i < dim; i++)
			{
				if (sigma_l[i] >= 0 && sigma_r[i] >= 0)
				{
					/* estimates are okay */
					double var = pow (sigma_l[i], beta) + pow (sigma_r[i], beta);
					
					if (var <= best_var)
					{
						found_best = 1;
						best_var = var;
						i_bisect = i;
						weight_l = pow (sigma_l[i], beta);
						weight_r = pow (sigma_r[i], beta);
					}
				}
				else
				{
					if (sigma_l[i] < 0)
					{
						GSL_ERROR ("no points in left-half space!", GSL_ESANITY);
					}
					if (sigma_r[i] < 0)
					{
						GSL_ERROR ("no points in right-half space!", GSL_ESANITY);
					}
				}
			}
		}
		
		if (!found_best)
		{
			/* All estimates were the same, so chose a direction at random */
			
			i_bisect = rng.RandUniformInt(dim);
		}
		
		xbi_l = xl[i_bisect];
		xbi_m = xmid[i_bisect];
		xbi_r = xu[i_bisect];
		
		/* Get the actual fractional sizes of the two "halves", and
		distribute the remaining calls among them */
		
		{
			double fraction_l = fabs ((xbi_m - xbi_l) / (xbi_r - xbi_l));
			double fraction_r = 1 - fraction_l;
			
			double a = fraction_l * weight_l;
			double b = fraction_r * weight_r;
			
			calls_l = min_calls + (calls - 2 * min_calls) * a / (a + b);
			calls_r = min_calls + (calls - 2 * min_calls) * b / (a + b);
		}
		
		/* Compute the integral for the left hand side of the bisection */
		
		/* Due to the recursive nature of the algorithm we must allocate
		some new memory for each recursive call */
		
		{
			
			double *xu_tmp = new double[dim];
			
			if (xu_tmp == 0)
			{
				GSL_ERROR_VAL ("out of memory for left workspace", GSL_ENOMEM, 0);
			}
			
			for (i = 0; i < dim; i++)
			{
				xu_tmp[i] = xu[i];
			}
			
			xu_tmp[i_bisect] = xbi_m;
			
			res_l  = MiserIntegrate(xl, xu_tmp, dim, calls_l,err_l);
			delete[] xu_tmp;
			xu_tmp = NULL;
			return res_l;
		}
		
		/* Compute the integral for the right hand side of the bisection */
		
		{
			
			double *xl_tmp = new double[dim];
			
			if (xl_tmp == 0)
			{
				GSL_ERROR_VAL ("out of memory for right workspace", GSL_ENOMEM, 0);
			}
			
			for (i = 0; i < dim; i++)
			{
				xl_tmp[i] = xl[i];
			}
			
			xl_tmp[i_bisect] = xbi_m;
			
			res_r = MiserIntegrate(xl_tmp, xu, dim, calls_r,err_r);

			delete [] xl_tmp;
			
			return res_r;
		}
		
		
		abserr = sqrt (err_l * err_l + err_r * err_r);
		
		return (res_l + res_r);
	}
	
	double CMonteCarlo::EstimateCorrmc( const double xl[], const double xu[], size_t dim, size_t calls, double &abserr, double sigma_l[], double sigma_r[] )
	{
		size_t i, n;
		monte_miser_state* state = (monte_miser_state*)mc_state;
		double *x = state->x;
		double *fsum_l = state->fsum_l;
		double *fsum_r = state->fsum_r;
		double *fsum2_l = state->fsum2_l;
		double *fsum2_r = state->fsum2_r;
		double *xmid = state->xmid;
		size_t *hits_l = state->hits_l;
		size_t *hits_r = state->hits_r;
		
		double m = 0.0, q = 0.0; 
		double vol = 1.0;
		
		for (i = 0; i < dim; i++)
		{
			vol *= xu[i] - xl[i];
			hits_l[i] = hits_r[i] = 0;
			fsum_l[i] = fsum_r[i] = 0.0;
			fsum2_l[i] = fsum2_r[i] = 0.0;
			sigma_l[i] = sigma_r[i] = -1;
		}
		
		for (n = 0; n < calls; n++)
		{
			double fval;
			
			unsigned int j = (n/2) % dim;
			unsigned int side = (n % 2);
			
			for (i = 0; i < dim; i++)
			{
				double z = rng.RandUniform(true);
				
				if (i != j) 
				{
					x[i] = xl[i] + z * (xu[i] - xl[i]);
				}
				else
				{
					if (side == 0) 
					{
						x[i] = xmid[i] + z * (xu[i] - xmid[i]);
					}
					else
					{
						x[i] = xl[i] + z * (xmid[i] - xl[i]);
					}
				}
			}
			
			fval = MONTE_FN_EVAL (mc_func, x);
			
			/* recurrence for mean and variance */
			{
				double d = fval - m;
				m += d / (n + 1.0);
				q += d * d * (n / (n + 1.0));
			}
			
			/* compute the variances on each side of the bisection */
			for (i = 0; i < dim; i++)
			{
				if (x[i] <= xmid[i])
				{
					fsum_l[i] += fval;
					fsum2_l[i] += fval * fval;
					hits_l[i]++;
				}
				else
				{
					fsum_r[i] += fval;
					fsum2_r[i] += fval * fval;
					hits_r[i]++;
				}
			}
		}
		
		for (i = 0; i < dim; i++)
		{
			double fraction_l = (xmid[i] - xl[i]) / (xu[i] - xl[i]);
			
			if (hits_l[i] > 0)
			{
				fsum_l[i] /= hits_l[i];
				sigma_l[i] = sqrt (fsum2_l[i] - fsum_l[i] * fsum_l[i] / hits_l[i]);
				sigma_l[i] *= fraction_l * vol / hits_l[i];
			}
			
			if (hits_r[i] > 0)
			{
				fsum_r[i] /= hits_r[i];
				sigma_r[i] = sqrt (fsum2_r[i] - fsum_r[i] * fsum_r[i] / hits_r[i]);
				sigma_r[i] *= (1 - fraction_l) * vol / hits_r[i];
			}
		}
		
		if (calls < 2)
		{
			abserr = CMath::POSINF;
		}
		else
		{
			abserr = vol * sqrt (q / (calls * (calls - 1.0)));
		}
		
		return vol * m;	
	}
	
	double CMonteCarlo::VegasIntegrate( double* xl, double* xu, const size_t dim, size_t calls, double& abserr )
	{
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		double cum_int, cum_sig;
		size_t i, k, it;
		
		if (dim != state->dim)
		{
			GSL_ERROR ("number of dimensions must match allocated size", GSL_EINVAL);
		}
		
		for (i = 0; i < dim; i++)
		{
			if (xu[i] <= xl[i])
			{
				GSL_ERROR ("xu must be greater than xl", GSL_EINVAL);
			}
			
			if (xu[i] - xl[i] > DOUBLE_MAX)
			{
				GSL_ERROR ("Range of integration is too large, please rescale",
					GSL_EINVAL);
			}
		}
		
		if (state->stage == 0)
		{
			InitGrid(xl, xu, dim);
			
			if (state->verbose >= 0)
			{
				PrintLimit(xl, xu, dim);
			}
		}
		
		if (state->stage <= 1)
		{
			state->wtd_int_sum = 0;
			state->sum_wgts = 0;
			state->chi_sum = 0;
			state->it_num = 1;
			state->samples = 0;
		}
		
		if (state->stage <= 2)
		{
			unsigned int bins = state->bins_max;
			unsigned int boxes = 1;
			
			if (state->mode != VEGAS_MODE_IMPORTANCE_ONLY)
			{
				/* shooting for 2 calls/box */
				
				boxes = floor (pow (calls / 2.0, 1.0 / dim));
				state->mode = VEGAS_MODE_IMPORTANCE;
				
				if (2 * boxes >= state->bins_max)
				{
					/* if bins/box < 2 */
					int box_per_bin = MAX(boxes / state->bins_max, 1);
					
					bins = MIN(boxes / box_per_bin, state->bins_max);
					boxes = box_per_bin * bins;
					
					state->mode = VEGAS_MODE_STRATIFIED;
				}
			}
			
			{
				double tot_boxes = pow ((double) boxes, (double) dim);
				state->calls_per_box = MAX(calls / tot_boxes, 2);
				calls = state->calls_per_box * tot_boxes;
			}
			
			/* total volume of x-space/(avg num of calls/bin) */
			state->jac = state->vol * pow ((double) bins, (double) dim) / calls;
			
			state->boxes = boxes;
			
			/* If the number of bins changes from the previous invocation, bins
			are expanded or contracted accordingly, while preserving bin
			density */
			
			if (bins != state->bins)
			{
				ResizeGrid(bins);
				
				if (state->verbose > 1)
				{
					PrintGrid(dim);
				}
			}
			
			if (state->verbose >= 0)
			{
				PrintHead(dim, calls, state->it_num, state->bins, state->boxes);
			}
		}
		
		state->it_start = state->it_num;
		
		cum_int = 0.0;
		cum_sig = 0.0;
		
		state->chisq = 0.0;
		
		for (it = 0; it < state->iterations; it++)
		{
			double intgrl = 0.0, intgrl_sq = 0.0;
			double sig = 0.0;
			double wgt;
			size_t calls_per_box = state->calls_per_box;
			double jacbin = state->jac;
			double *x = state->x;
			int *bin = state->bin;
			
			state->it_num = state->it_start + it;
			
			ResetGridValues();
			InitBoxCoord(state->box);
			
			do
			{
				double m = 0, q = 0;
				double f_sq_sum = 0.0;
				
				for (k = 0; k < calls_per_box; k++)
				{
					double fval, bin_vol;
					
					RandomPoint(x, bin, &bin_vol, state->box, xl, xu);
					
					fval = jacbin * bin_vol * MONTE_FN_EVAL (mc_func, x);
					
					/* recurrence for mean and variance */
					
					{
						double d = fval - m;
						m += d / (k + 1.0);
						q += d * d * (k / (k + 1.0));
					}
					
					if (state->mode != VEGAS_MODE_STRATIFIED)
					{
						double f_sq = fval * fval;
						AccumulateDistribution(bin, f_sq);
					}
				}
				
				intgrl += m * calls_per_box;
				
				f_sq_sum = q * calls_per_box ;
				
				sig += f_sq_sum ;
				
				if (state->mode == VEGAS_MODE_STRATIFIED)
				{
					AccumulateDistribution(bin, f_sq_sum);
				}
			}
			while (ChangeBoxCoord(state->box));
			
			/* Compute final results for this iteration   */
			
			sig = sig / (calls_per_box - 1.0)  ;
			
			if (sig > 0) 
			{
				wgt = 1.0 / sig;
			}
			else if (state->sum_wgts > 0) 
			{
				wgt = state->sum_wgts / state->samples;
			}
			else 
			{
				wgt = 0.0;
			}
			
			intgrl_sq = intgrl * intgrl;
			
			state->result = intgrl;
			state->sigma  = sqrt(sig);
			
			if (wgt > 0.0)
			{
				state->samples++ ;
				state->sum_wgts += wgt;
				state->wtd_int_sum += intgrl * wgt;
				state->chi_sum += intgrl_sq * wgt;
				
				cum_int = state->wtd_int_sum / state->sum_wgts;
				cum_sig = sqrt (1 / state->sum_wgts);
				
				if (state->samples > 1)
				{
					state->chisq = (state->chi_sum - state->wtd_int_sum * cum_int) /
						(state->samples - 1.0);
				}
			}
			else
			{
				cum_int += (intgrl - cum_int) / (it + 1.0);
				cum_sig = 0.0;
			}         
			
			
			if (state->verbose >= 0)
			{
				PrintRes
					(
					state->it_num, intgrl, sqrt (sig), cum_int, cum_sig,
					state->chisq
					);
				if (it + 1 == state->iterations && state->verbose > 0)
				{
					PrintGrid(dim);
				}
			}
			
			if (state->verbose > 1)
			{
				PrintDist(dim);
			}
			
			RefineGrid();
			
			if (state->verbose > 1)
			{
				PrintGrid(dim);
			}
			
		}
	
		/* By setting stage to 1 further calls will generate independent
		estimates based on the same grid, although it may be rebinned. */
		
		state->stage = 1;  
		
		
		abserr = cum_sig;
		
		return cum_int;
	}
	
	
	void CMonteCarlo::PrintLimit( double xl[], double xu[], unsigned long dim )
	{
		unsigned long j;
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		fprintf (state->ostream, "The limits of integration are:\n");
		for (j = 0; j < dim; ++j)
			fprintf (state->ostream, "\nxl[%lu]=%f    xu[%lu]=%f", j, xl[j], j, xu[j]);
		fprintf (state->ostream, "\n");
		fflush (state->ostream);		
	}
	
	void CMonteCarlo::PrintHead( unsigned long num_dim, unsigned long calls, unsigned int it_num, unsigned int bins, unsigned int boxes )
	{
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		fprintf (state->ostream,
			"\nnum_dim=%lu, calls=%lu, it_num=%d, max_it_num=%d ",
			num_dim, calls, it_num, state->iterations);
		fprintf (state->ostream,
			"verb=%d, alph=%.2f,\nmode=%d, bins=%d, boxes=%d\n",
			state->verbose, state->alpha, state->mode, bins, boxes);
		fprintf (state->ostream,
			"\n       single.......iteration                   ");
		fprintf (state->ostream, "accumulated......results   \n");
		
		fprintf (state->ostream,
			"iteration     integral    sigma             integral   ");
		fprintf (state->ostream, "      sigma     chi-sq/it\n\n");
		fflush (state->ostream);		
	}
	
	void CMonteCarlo::PrintRes( unsigned int itr, double res, double err, double cum_res, double cum_err, double chi_sq )
	{
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		fprintf (state->ostream,
			"%4d        %6.4e %10.2e          %6.4e      %8.2e  %10.2e\n",
			itr, res, err, cum_res, cum_err, chi_sq);
		fflush (state->ostream);		
	}
	
	void CMonteCarlo::PrintDist( unsigned long dim )
	{
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		unsigned long i, j;
		int p = state->verbose;
		if (p < 1)
			return;
		
		for (j = 0; j < dim; ++j)
		{
			fprintf (state->ostream, "\n axis %lu \n", j);
			fprintf (state->ostream, "      x   g\n");
			for (i = 0; i < state->bins; i++)
			{
				fprintf (state->ostream, "weight [%11.2e , %11.2e] = ", 
					COORD (state, i, j), COORD(state,i+1,j));
				fprintf (state->ostream, " %11.2e\n", VALUE (state, i, j));
				
			}
			fprintf (state->ostream, "\n");
		}
		fprintf (state->ostream, "\n");
		fflush (state->ostream);	
	}
	
	void CMonteCarlo::PrintGrid( unsigned long dim )
	{
		monte_vegas_state* state = (monte_vegas_state*)mc_state;
		unsigned long i, j;
		int p = state->verbose;
		if (p < 1)
			return;
		
		for (j = 0; j < dim; ++j)
		{
			fprintf (state->ostream, "\n axis %lu \n", j);
			fprintf (state->ostream, "      x   \n");
			for (i = 0; i <= state->bins; i++)
			{
				fprintf (state->ostream, "%11.2e", COORD (state, i, j));
				if (i % 5 == 4)
					fprintf (state->ostream, "\n");
			}
			fprintf (state->ostream, "\n");
		}
		fprintf (state->ostream, "\n");
		fflush (state->ostream);	
	}
	
	void CMonteCarlo::RefineGrid()
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;
		size_t i, j, k;
		size_t dim = s->dim;
		size_t bins = s->bins;
		
		for (j = 0; j < dim; j++)
		{
			double grid_tot_j, tot_weight;
			double * weight = s->weight;
			
			double oldg = VALUE (s, 0, j);
			double newg = VALUE (s, 1, j);
			
			VALUE (s, 0, j) = (oldg + newg) / 2;
			grid_tot_j = VALUE (s, 0, j);
			
			/* This implements gs[i][j] = (gs[i-1][j]+gs[i][j]+gs[i+1][j])/3 */
			
			for (i = 1; i < bins - 1; i++)
			{
				double rc = oldg + newg;
				oldg = newg;
				newg = VALUE (s, i + 1, j);
				VALUE (s, i, j) = (rc + newg) / 3;
				grid_tot_j += VALUE (s, i, j);
			}
			VALUE (s, bins - 1, j) = (newg + oldg) / 2;
			
			grid_tot_j += VALUE (s, bins - 1, j);
			
			tot_weight = 0;
			
			for (i = 0; i < bins; i++)
			{
				weight[i] = 0;
				
				if (VALUE (s, i, j) > 0)
				{
					oldg = grid_tot_j / VALUE (s, i, j);
					/* damped change */
					weight[i] = pow (((oldg - 1) / oldg / log (oldg)), s->alpha);
				}
				
				tot_weight += weight[i];
				
#ifdef DEBUG
				printf("weight[%d] = %g\n", i, weight[i]);
#endif
			}
			
			{
				double pts_per_bin = tot_weight / bins;
				
				double xold;
				double xnew = 0;
				double dw = 0;
				i = 1;
				
				for (k = 0; k < bins; k++)
				{
					dw += weight[k];
					xold = xnew;
					xnew = COORD (s, k + 1, j);
					
					for (; dw > pts_per_bin; i++)
					{
						dw -= pts_per_bin;
						NEW_COORD (s, i) = xnew - (xnew - xold) * dw / weight[k];
					}
				}
				
				for (k = 1 ; k < bins ; k++)
				{
					COORD(s, k, j) = NEW_COORD(s, k);
				}
				
				COORD (s, bins, j) = 1;
			}
		}	
	}
	
	void CMonteCarlo::ResizeGrid( unsigned int bins )
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;
		size_t j, k;
		size_t dim = s->dim;
		
		/* weight is ratio of bin sizes */
		
		double pts_per_bin = (double) s->bins / (double) bins;
		
		for (j = 0; j < dim; j++)
		{
			double xold;
			double xnew = 0;
			double dw = 0;
			int i = 1;
			
			for (k = 1; k <= s->bins; k++)
			{
				dw += 1.0;
				xold = xnew;
				xnew = COORD (s, k, j);
				
				for (; dw > pts_per_bin; i++)
				{
					dw -= pts_per_bin;
					NEW_COORD (s, i) = xnew - (xnew - xold) * dw;
				}
			}
			
			for (k = 1 ; k < bins; k++)
			{
				COORD(s, k, j) = NEW_COORD(s, k);
			}
			
			COORD (s, bins, j) = 1;
		}
		
		s->bins = bins;	
	}
	
	void CMonteCarlo::RandomPoint( double x[], int bin[], double *bin_vol, const int box[], const double xl[], const double xu[] )
	{
		
	/* Use the random number generator r to return a random position x
	in a given box.  The value of bin gives the bin location of the
		random position (there may be several bins within a given box) */
		monte_vegas_state* s = (monte_vegas_state*)mc_state;
		double vol = 1.0;
		
		size_t j;
		
		size_t dim = s->dim;
		size_t bins = s->bins;
		size_t boxes = s->boxes;
		
		//DISCARD_POINTER(xu); /* prevent warning about unused parameter */
		
		for (j = 0; j < dim; ++j)
		{
		/* box[j] + ran gives the position in the box units, while z
			is the position in bin units.  */
			
			double z = ((box[j] + rng.RandUniform(true)) / boxes) * bins;
			
			int k = z;
			
			double y, bin_width;
			
			bin[j] = k;
			
			if (k == 0)
			{
				bin_width = COORD (s, 1, j);
				y = z * bin_width;
			}
			else
			{
				bin_width = COORD (s, k + 1, j) - COORD (s, k, j);
				y = COORD (s, k, j) + (z - k) * bin_width;
			}
			
			x[j] = xl[j] + y * s->delx[j];
			
			vol *= bin_width;
		}
		
		*bin_vol = vol;		
	}
	
	void CMonteCarlo::AccumulateDistribution( int bin[], double y )
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;
		size_t j;
		size_t dim = s->dim;
		
		for (j = 0; j < dim; j++)
		{
			int i = bin[j];
			VALUE (s, i, j) += y;
		}		
	}
	
	void CMonteCarlo::ResetGridValues()
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;
		size_t i, j;
		
		size_t dim = s->dim;
		size_t bins = s->bins;
		
		for (i = 0; i < bins; i++)
		{
			for (j = 0; j < dim; j++)
			{
				VALUE (s, i, j) = 0.0;
			}
		}	
	}
	
	void CMonteCarlo::InitGrid( const double xl[], const double xu[], size_t dim )
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;

		size_t j;
		
		double vol = 1.0;
		
		s->bins = 1;
		
		for (j = 0; j < dim; j++)
		{
			double dx = xu[j] - xl[j];
			s->delx[j] = dx;
			vol *= dx;
			
			COORD (s, 0, j) = 0.0;
			COORD (s, 1, j) = 1.0;
		}
		
		s->vol = vol;		
	}
	
	int CMonteCarlo::ChangeBoxCoord( int box[] )
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;

		int j = s->dim - 1;
		
		int ng = s->boxes;
		
		while (j >= 0)
		{
			box[j] = (box[j] + 1) % ng;
			
			if (box[j] != 0)
			{
				return 1;
			}
			
			j--;
		}
		
		return 0;	
	}
	
	void CMonteCarlo::InitBoxCoord( int box[] )
	{
		monte_vegas_state* s = (monte_vegas_state*)mc_state;

		size_t i;
		
		size_t dim = s->dim;
		
		for (i = 0; i < dim; i++)
		{
			box[i] = 0;
		}	
	}
	
	inline void* CMonteCarlo::GetState() const
	{
		return mc_state;	
	}
	
#undef  GSL_MONTE_FN_EVAL
}