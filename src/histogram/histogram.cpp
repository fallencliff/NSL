/********************************************************************
	filename: 	histogram.cpp
	author:		huzhijian
	created:	5:5:2010   16:43
	brief:	Öù×´Í¼
*********************************************************************/

#include <math.h>
#include <gsl_errno.h> 
#include <histogram.h>

namespace gslcpp
{
	CHistogram::CHistogram(const size_t n, bool b_zero)
		:gsldata(NULL), gslpdf(NULL)
	{
		gsldata = (b_zero == false) ? Alloc(n) : Calloc(n);
	}
	
	CHistogram::CHistogram( const size_t n, const double xmin, const double xmax )
		:gsldata(NULL), gslpdf(NULL)
	{
		
		if (xmin >= xmax)
		{
			gsl_error("xmin must be less than xmax", __FILE__, __LINE__, GSL_EINVAL);
		}
		
		gsldata = Calloc (n);
		
		if (gsldata == 0)
		{
			gsl_error("memory alloc failed ", __FILE__, __LINE__, GSL_FAILURE);
		}
		
		make_uniform (n, xmin, xmax);
		
		return;	
	}
	
	gsl_histogram* CHistogram::Alloc(const size_t n)
	{
		gsl_histogram *h;
		
		if (n == 0)
		{
			gsl_error("histogram length n must be positive integer", __FILE__, __LINE__,
				GSL_EDOM);
		}
		
		h = (gsl_histogram *) new gsl_histogram;
		
		if (h == 0)
		{
			gsl_error ("failed to allocate space for histogram struct",__FILE__, __LINE__,
				GSL_ENOMEM);
		}
		
		h->range = (double *) new double[n+1];
		
		if (h->range == 0)
		{
			delete h;         /* exception in constructor, avoid memory leak */
			
			gsl_error ("failed to allocate space for histogram ranges", __FILE__, __LINE__,
				GSL_ENOMEM);
		}
		
		h->bin = (double *) new double[n];
		
		if (h->bin == 0)
		{
			delete [] h->range;
			delete h;         /* exception in constructor, avoid memory leak */
			
			gsl_error ("failed to allocate space for histogram bins", __FILE__, __LINE__,
				GSL_ENOMEM);
		}
		
		h->n = n;

		return h;
	}
	
	void CHistogram::Free()
	{
		if (gsldata)
		{
			if (gsldata->range)
			{
				delete [] gsldata->range;
				gsldata->range = NULL;
			}
			if (gsldata->bin)
			{
				delete[] gsldata->bin;
				gsldata->bin = NULL;
			}
			delete gsldata;
			gsldata = NULL;
		}

		if (gslpdf)
		{
			if (gslpdf->range)
			{
				delete [] gslpdf->range;
			}
			if (gslpdf->sum)
			{
				delete [] gslpdf->sum;
			}
			delete gslpdf;
			gslpdf = NULL;
		}
			
	}
	
	int CHistogram::SetRanges( const double range[], const size_t size )
	{
		size_t i;
		const size_t n = gsldata->n;
		
		if (size != (n+1))
		{
			GSL_ERROR ("size of range must match size of histogram", GSL_EINVAL);
		}
		
		/* initialize ranges */
		
		for (i = 0; i <= n; i++)
		{
			gsldata->range[i] = range[i];
		}
		
		/* clear contents */
		
		for (i = 0; i < n; i++)
		{
			gsldata->bin[i] = 0;
		}
		
		return GSL_SUCCESS;		
	}
	
	int CHistogram::SetRanges( const double xmin, const double xmax )
	{
		size_t i;
		const size_t n = gsldata->n;
		
		if (xmin >= xmax)
		{
			GSL_ERROR ("xmin must be less than xmax", GSL_EINVAL);
		}
		
		/* initialize ranges */
		
		make_uniform (n, xmin, xmax);
		
		/* clear contents */
		
		for (i = 0; i < n; i++)
		{
			gsldata->bin[i] = 0;
		}
		
		return GSL_SUCCESS;	
	}
	
	void CHistogram::make_uniform(const size_t n,const double xmin, const double xmax)
	{
		size_t i;
		
		for (i = 0; i <= n; i++)
		{
			double f1 = ((double) (n-i) / (double) n);
			double f2 = ((double) i / (double) n);
			gsldata->range[i] = f1 * xmin +  f2 * xmax;
		}		
	}
	
	gsl_histogram* CHistogram::Calloc( const size_t n )
	{
		gsl_histogram * h = Alloc (n);
		
		if (h == 0)
		{
			return h;
		}
		
		{
			size_t i;
			
			for (i = 0; i < n + 1; i++)
			{
				h->range[i] = i;
			}
			
			for (i = 0; i < n; i++)
			{
				h->bin[i] = 0;
			}
		}
		
		h->n = n;
		
		return h;		
	}
	
	void CHistogram::Resize( const size_t n )
	{
		if (n == gsldata->n)
		{
			return;
		}
		Free();
		gsldata = Alloc(n);
	}
	
	CHistogram& CHistogram::operator=( const CHistogram & src )
	{
		size_t n = src.gsldata->n;

		Resize(n);
		size_t i;
		for (i = 0; i <= n; i++)
		{
			gsldata->range[i] = src.gsldata->range[i];
		}
		
		for (i = 0; i < n; i++)
		{
			gsldata->bin[i] = src.gsldata->bin[i];
		}
		
		return *this;		
	}
	
	double CHistogram::Get( const size_t index ) const
	{
	
		if (index >= gsldata->n)
		{
			GSL_ERROR_VAL ("index lies outside valid range of 0 .. n - 1",
				GSL_EDOM, 0);
		}
		
		return gsldata->bin[index];	
	}
	
	int CHistogram::GetRange( const size_t index, double &lower, double&upper ) const
	{
		if (index >= gsldata->n)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. n - 1", GSL_EDOM);
		}
		
		lower = gsldata->range[index];
		upper = gsldata->range[index + 1];
		
		return GSL_SUCCESS;		
	}
	
	size_t CHistogram::Find(const size_t n, const double range[], const double x ) const
	{
		size_t lower, upper, mid, res;
		//const double* range = gsldata->range;;
		if (x < range[0] || x >= range[n])
		{
			GSL_ERROR ("x not found in range of h", GSL_EDOM);
		}

		/* optimize for linear case */
		
#ifdef LINEAR_OPT
		{
			double u =  (x - range[0]) / (range[n] - range[0]);
			res = (size_t) (u * n);
		}
		
		if (x >= range[res] && x < range[res + 1])
		{
			return res;
		}
#endif
		
		/* perform binary search */
		
		upper = n ;
		lower = 0 ;
		
		while (upper - lower > 1)
		{
			mid = (upper + lower) / 2 ;
			
			if (x >= range[mid])
			{
				lower = mid ;
			}
			else
			{
				upper = mid ;
			}
		}
		
		res = lower ;
		
		/* sanity check the result */
		
		if (x < range[lower] || x >= range[lower + 1])
		{
			GSL_ERROR ("x not found in range", GSL_ESANITY);
		}

		return res;
	
	}
	inline double CHistogram::Max() const
	{
		return gsldata->range[gsldata->n];	
	}

	inline double CHistogram::Min() const
	{
		return gsldata->range[0];
	}
	
	inline double CHistogram::Bins() const
	{
		return gsldata->n;		
	}
	
	void CHistogram::Reset()
	{
		const size_t n = gsldata->n;
		double* bin = gsldata->bin;
		for (size_t i = 0; i < n; i++)
		{
			bin[i] = 0;
		}
	}
	
	int CHistogram::Increment( const double x )
	{
		return Increment (x, 1.0);	
	}
	
	int CHistogram::Increment( const double x, const double weight )
	{
		size_t index = Find(gsldata->n, gsldata->range, x);
		
		if (index >= gsldata->n)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. n - 1",
				GSL_ESANITY);
		}
		
		gsldata->bin[index] += weight;
		
		return GSL_SUCCESS;		
	}
	
	double CHistogram::MaxBin( size_t *index /*= NULL*/ ) const
	{
		size_t imax = 0;
		const double* bin = gsldata->bin;
		const size_t n = gsldata->n;
		double max = bin[0];

		for (size_t i = 0; i < n; i++)
		{
			if (bin[i] > max)
			{
				max = bin[i];
				imax = i;
			}
		}
		if (index != NULL)
		{
			*index = imax;
		}
		return max;		
	}
	
	double CHistogram::MinBin( size_t *index /*= NULL*/ ) const
	{
		size_t i;
		size_t imin = 0;
		const double* bin = gsldata->bin;
		const size_t n = gsldata->n;
		double min = bin[0];
		for (i = 0; i < n; i++)
		{
			if (bin[i] < min)
			{
				min = bin[i];
				imin = i;
			}
		}
		if (index != NULL)
		{
			*index = imin;
		}
		return min;		
	}
	
	double CHistogram::Sum() const
	{
		double sum = 0;
		size_t i=0;
		size_t n = gsldata->n;
		const double* bin = gsldata->bin;

		while(i < n)
			sum += bin[i++];
		
		return sum;		
	}
	
	double CHistogram::Mean() const
	{
		const size_t n = gsldata->n;
		
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wmean = 0;
		long double W = 0;
		const double *range = gsldata->range;
		const double *bin = gsldata->bin;

		for (size_t i = 0; i < n; i++)
		{
			double xi = (range[i + 1] + range[i]) / 2;
			double wi = bin[i];
			
			if (wi > 0)
			{
				W += wi;
				wmean += (xi - wmean) * (wi / W);
			}
		}
		
		return wmean;		
	}
	
	double CHistogram::Sigma() const
	{
		const size_t n = gsldata->n;
		size_t i;
		
		long double wvariance = 0 ;
		long double wmean = 0;
		long double W = 0;
		const double *range = gsldata->range;
		const double *bin = gsldata->bin;
		/* FIXME: should use a single pass formula here, as given in
		N.J.Higham 'Accuracy and Stability of Numerical Methods', p.12 */
		
		/* Compute the mean */
		
		for (i = 0; i < n; i++)
		{
			double xi = (range[i + 1] + range[i]) / 2;
			double wi = bin[i];
			
			if (wi > 0)
			{
				W += wi;
				wmean += (xi - wmean) * (wi / W);
			}
		}
		
		/* Compute the variance */
		
		W = 0.0;
		
		for (i = 0; i < n; i++)
		{
			double xi = ((range[i + 1]) + (range[i])) / 2;
			double wi = bin[i];
			
			if (wi > 0) {
				const long double delta = (xi - wmean);
				W += wi ;
				wvariance += (delta * delta - wvariance) * (wi / W);
			}
		}
		
		{
			double sigma = sqrt(wvariance);
			return sigma;
		}		
	}
	
	bool CHistogram::operator==( const CHistogram& other ) const
	{
		if (gsldata->n != other.gsldata->n)
		{
			return false;
		}
		
		{
			size_t i;
			const size_t n = gsldata->n;
			/* init ranges */
			const double *range = gsldata->range;
			const double *range_other = other.gsldata->range;
			for (i = 0; i <= n; i++)
			{
				if (range[i] != range_other[i])
				{
					return false;
				}
			}
		}
		
		return true;		
	}
	
	bool CHistogram::operator!=( const CHistogram& other ) const
	{
		return !(*this == other);		
	}
	
	CHistogram& CHistogram::operator+=( const CHistogram & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t n  = gsldata->n;

		for (size_t i = 0; i < n; i++)
		{
			bin[i] += bin_other[i];
		}
		
		return *this;
	}
	
	CHistogram& CHistogram::operator+=( const double & other )
	{
		double *bin = gsldata->bin;
		const size_t n  = gsldata->n;
		
		for (size_t i = 0; i < n; i++)
		{
			bin[i] += other;
		}
		
		return *this;	
	}
	
	CHistogram& CHistogram::operator-=( const CHistogram & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t n  = gsldata->n;
		
		for (size_t i = 0; i < n; i++)
		{
			bin[i] -= bin_other[i];
		}
		
		return *this;	
	}
	
	CHistogram& CHistogram::operator*=( const CHistogram & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t n  = gsldata->n;
		
		for (size_t i = 0; i < n; i++)
		{
			bin[i] *= bin_other[i];
		}
		
		return *this;	
	}
	
	CHistogram& CHistogram::operator*=( const double & other )
	{
		double *bin = gsldata->bin;
		const size_t n  = gsldata->n;
		
		for (size_t i = 0; i < n; i++)
		{
			bin[i] *= other;
		}
		
		return *this;
	}
	
	CHistogram& CHistogram::operator/=( const CHistogram & other )
	{
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t n  = gsldata->n;
		
		for (size_t i = 0; i < n; i++)
		{
			bin[i] /= bin_other[i];
		}
		
		return *this;		
	}
	
	int CHistogram::Save( const char* filename ) const
	{
		const size_t n = gsldata->n;
		const double* bin = gsldata->bin;
		const double* range = gsldata->range;
		FILE* fp = fopen(filename, "w");

		for (size_t i =0; i<n; i++)
		{
			if (fprintf(fp, "%g\t%g\t%g\n", range[i], range[i+1], bin[i]) < 0)
			{
				fclose(fp);
				GSL_ERROR ("Histogram save failed", GSL_EFAILED);
			}
		}
		fclose(fp);
		return GSL_SUCCESS;		
	}
	
	int CHistogram::Load( const char* filename )
	{
		const size_t n = gsldata->n;
		double* bin = gsldata->bin;
		double* range = gsldata->range;
		double upper;
		FILE* fp = fopen(filename, "r");
		
		for (size_t i =0; i<n; i++)
		{
			if (fscanf(fp, "%lg\t%lg\t%lg\n", range+i, &upper, bin+i) != 3)
			{
				fclose(fp);
				GSL_ERROR ("Histogram load failed", GSL_EFAILED);
			}
		}
		fclose(fp);
		range[n] = upper;

		return GSL_SUCCESS;
	}
	
	gsl_histogram_pdf* CHistogram::AllocPdf( const size_t n )
	{
		gsl_histogram_pdf *p;
		
		if (n == 0)
		{
			GSL_ERROR_VAL ("histogram pdf length n must be positive integer",
				GSL_EDOM, 0);
		}
		
		p = (gsl_histogram_pdf *) new gsl_histogram_pdf;
		
		if (p == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for histogram pdf struct",
				GSL_ENOMEM, 0);
		}
		
		p->range = (double *) new double[n+1];
		
		if (p->range == 0)
		{
			delete p;         /* exception in constructor, avoid memory leak */
			p = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram pdf ranges",
				GSL_ENOMEM, 0);
		}
		
		p->sum = (double *) new double[n+1];
		
		if (p->sum == 0)
		{
			delete [] p->range;
			delete p;         /* exception in constructor, avoid memory leak */
			p = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram pdf sums",
				GSL_ENOMEM, 0);
		}
		
		p->n = n;
		
		return p;		
	}
	
	int CHistogram::InitPdf()
	{
		size_t i;
		size_t n = gslpdf->n;
		const double* bin = gsldata->bin;
		const double* range = gsldata->range;
		double* range_pdf = gslpdf->range;
		if (n != gsldata->n)
		{
			GSL_ERROR ("histogram length must match pdf length", GSL_EINVAL);
		}
		
		for (i = 0; i < n; i++)
		{
			if (bin[i] < 0)
			{
				GSL_ERROR ("histogram bins must be non-negative to compute"
					"a probability distribution", GSL_EDOM);
			}
		}
		
		for (i = 0; i < n + 1; i++)
		{
			range_pdf[i] = range[i];
		}
		
		{
			double mean = 0, sum = 0;
			
			for (i = 0; i < n; i++)
			{
				mean += (bin[i] - mean) / ((double) (i + 1));
			}
			
			gslpdf->sum[0] = 0;
			
			for (i = 0; i < n; i++)
			{
				sum += (bin[i] / mean) / n;
				gslpdf->sum[i + 1] = sum;
			}
		}
		
		return GSL_SUCCESS;		
	}
	
	double CHistogram::PdfSample( double r )
	{
		if (!gslpdf)
		{
			gslpdf = AllocPdf(gsldata->n);
			InitPdf();
		}
		
		/* Wrap the exclusive top of the bin down to the inclusive bottom of
		the bin. Since this is a single point it should not affect the
		distribution. */
		
		if (r == 1.0)
		{
			r = 0.0;
		}
		
		size_t i = Find (gslpdf->n, gslpdf->sum, r);
		
		
		double delta = (r - gslpdf->sum[i]) / (gslpdf->sum[i + 1] - gslpdf->sum[i]);
		double x = gslpdf->range[i] + delta * (gslpdf->range[i + 1] - gslpdf->range[i]);
		return x;
		
	}
	
	CHistogram::~CHistogram()
	{
		Free();		
	}
}