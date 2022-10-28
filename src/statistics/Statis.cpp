/********************************************************************
	filename: 	Statis.cpp
	author:		huzhijian
	created:	5:5:2010   16:43
	brief:	
*********************************************************************/

#include <Statis.h>
#include <MathConst.h>
#include <math.h>


namespace gslcpp
{
	CStatis::CStatis(const double * base, size_t n, size_t m_stride)
		: data(base), size(n), stride(m_stride), wdata(0), wsize(0), wstride(0)
	{}

	CStatis::~CStatis()
	{

	}

	double CStatis::get_mean() const
	{
	/* Compute the arithmetic mean of a dataset using the recurrence relation 
		mean_(n) = mean(n-1) + (data[n] - mean(n-1))/(n+1)   */
		
		double mean = 0;
		for (size_t i = 0; i < size; i++)
		{
			mean += (data[i * stride] - mean) / (i + 1);
		}
		return mean;
	}
	double CStatis::compute_variance(double mean) const
	{
		long double variance = 0 ;
		
		size_t i;

		/* find the sum of the squares */
		for (i = 0; i < size; i++)
		{
			const long double delta = (data[i * stride] - mean);
			variance += (delta * delta - variance) / (i + 1);
		}
		
		return variance;
	}

	double CStatis::variance() const
	{
		const double mean = get_mean();

		return variance(mean);
	}

	double CStatis::variance(double mean) const
	{
		double n = size;

		const double variance =  compute_variance(mean);

		return variance * ((double)n / (double)(n - 1));
	}

	double CStatis::std_deviation() const
	{
		const double mean = get_mean();

		return std_deviation(mean);
	}

	double CStatis::std_deviation(const double& mean) const
	{
		double n = size;

		const double variance = compute_variance(mean);

		const double sd = sqrt (variance * ((double)n / (double)(n - 1)));

		return sd;
	}
	
	double CStatis::douvariance_with_fixed_mean( const double & mean ) const
	{
		const double variance = compute_variance(mean);
		return variance;
	}
	
	double CStatis::std_deviation_with_fixed_mean( const double& mean ) const
	{
		const double variance = compute_variance(mean);
		return sqrt(variance);
	}
	
	double CStatis::absdev() const
{
		const double mean = get_mean();
		return absdev(mean);
	}
	
	double CStatis::absdev( const double& mean ) const
	{
		double sum = 0, absdev;
		size_t i;
		size_t n = size;
		/* find the sum of the absolute deviations */
		for (i = 0; i < n; i++)
		{
			const double delta = fabs(data[i * stride] - mean);
			sum += delta;
		}
		
		absdev = sum / n;
		
		return absdev;		
	}
	
	double CStatis::skew( const double& mean, const double& sd ) const
	{
		long double skew = 0;
		size_t i;

		/* find the sum of the cubed deviations, normalized by the sd. */
		
		/* we use a recurrence relation to stably update a running value so
		there aren't any large sums that can overflow */
		
		for (i = 0; i < size; i++)
		{
			const long double x = (data[i * stride] - mean) / sd;
			skew += (x * x * x - skew) / (i + 1);
		}
		
		return skew;		
	}
	
	double CStatis::skew() const
	{
		const double mean = get_mean();
		const double est_sd = std_deviation(mean);
		return skew(mean, est_sd);
	}
	
	double CStatis::kurtosis( const double& mean, const double& sd ) const
	{
		long double avg = 0, kurtosis;
		size_t i;
		/* find the fourth moment the deviations, normalized by the sd */
		
		/* we use a recurrence relation to stably update a running value so
		there aren't any large sums that can overflow */
		
		for (i = 0; i < size; i++)
		{
			const long double x = (data[i * stride] - mean) / sd;
			avg += (x * x * x * x - avg)/(i + 1);
		}
		
		kurtosis = avg - 3.0;  /* makes kurtosis zero for a Gaussian */
		
		return kurtosis;		
	}
	
	double CStatis::kurtosis() const
	{
		const double mean = get_mean();
		const double est_sd = std_deviation(mean);
		return kurtosis(mean, est_sd);
	}
	
	double CStatis::lag1_autocorrelation( const double& mean ) const
	{
		size_t i;
		long double r1 ;
		long double q = 0 ;
		long double v = (data[0 * stride] - mean) * (data[0 * stride] - mean) ;
		
		for (i = 1; i < size ; i++)
		{
			const long double delta0 = (data[(i-1) * stride] - mean);
			const long double delta1 = (data[i * stride] - mean);
			q += (delta0 * delta1 - q)/(i + 1);
			v += (delta1 * delta1 - v)/(i + 1);
		}
		
		r1 = q / v ;
		
		return r1;
	}
	
	double CStatis::lag1_autocorrelation() const
	{
		const double mean = get_mean();		
		return lag1_autocorrelation(mean);
	}
	
	double CStatis::covariance( const CStatis & other, const double& mean1, const double& mean2 ) const
	{
		size_t n = size;
		const double covariance = compute_covariance(other, mean1, mean2);
		return covariance * ((double)n / (double)(n - 1));
	}
	
	double CStatis::covariance( const CStatis& other ) const
	{
		const double mean1 = get_mean();
		const double mean2 = other.get_mean();
		return covariance(other, mean1, mean2);
	}
	
	double CStatis::compute_covariance( const CStatis& other , const double& mean1, const double& mean2) const
	{
		long double covariance = 0 ;
		
		size_t i;
		const double* data1 = data;
		const double* data2 = other.data;
		const size_t stride1 = stride;
		const size_t stride2 = other.stride;
		size_t n = size;
		/* find the sum of the squares */
		for (i = 0; i < n; i++)
		{
			const long double delta1 = (data1[i * stride1] - mean1);
			const long double delta2 = (data2[i * stride2] - mean2);
			covariance += (delta1 * delta2 - covariance) / (i + 1);
		}
		
		return covariance ;
	}
	
	void CStatis::AddWight( const double * base, size_t n, size_t stride /*= 1*/ )
	{
		wdata = const_cast<double*>(base);
		wsize = n;
		wstride = stride;
	}
	
	double CStatis::get_wmean() const
	{
	/* Compute the weighted arithmetic mean M of a dataset using the
	recurrence relation
	
	  M(n) = M(n-1) + (data[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
	  W(n) = W(n-1) + w(n)
	  
		*/
		
		double wmean = 0;
		double W = 0;
		
		size_t i;
		
		for (i = 0; i < size; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0)
			{
				W += wi;
				wmean += (data[i * stride] - wmean) * (wi / W);
			}
		}
		
		return wmean;	
	}
	
	double CStatis::wvariance() const
	{
		const double wmean = get_wmean();
		return wvariance(wmean);
	}	
	
	double CStatis::wvariance( const double& wmean ) const
	{
		const double variance = compute_wvariance(wmean);
		const double scale = compute_factor();
		
		return scale * variance;	
	}
	
	double CStatis::compute_wvariance( const double& wmean ) const
	{
		double wvariance = 0 ;
		double W = 0;
		
		size_t i;
		size_t n = size;
		/* find the sum of the squares */
		for (i = 0; i < n; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0) {
				const double delta = (data[i * stride] - wmean);
				W += wi ;
				wvariance += (delta * delta - wvariance) * (wi / W);
			}
		}
		
		return wvariance ;		
	}
	
	double CStatis::compute_factor() const
	{
		/* Find the factor ``N/(N-1)'' which multiplies the raw std dev */
		
		double a = 0 ;
		double b = 0;
		double factor;
		
		size_t i;
		
		/* find the sum of the squares */
		for (i = 0; i < size; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0)
			{
				a += wi ;
				b += wi * wi ;
			}
		}
		
		factor = (a*a) / ((a*a) - b);
		
		return factor ;	
	}
	
	double CStatis::wsd( const double& wmean ) const
	{
		const double variance = compute_wvariance(wmean);
		const double scale = compute_factor();
		const double wsd_ret = sqrt(scale * variance) ;
		
		return wsd_ret;		
	}
	
	double CStatis::wsd() const
	{
		const double wmean = get_wmean();
		return wsd(wmean);
	}
	
	double CStatis::wvariance_with_fixed_mean(const double& mean) const
	{
		const double wvariance = compute_wvariance(mean);		
		return wvariance;
	}
	
	double CStatis::wsd_with_fixed_mean( const double& mean ) const
	{
		const double wvariance	= compute_wvariance(mean);	
		return sqrt (wvariance);
	}
	
	double CStatis::wabsdev( const double& wmean ) const
	{
		double wabsdev = 0;
		double W = 0;
		
		size_t i;
		
		/* find the sum of the absolute deviations */
		for (i = 0; i < size; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0) {
				const double delta = fabs(data[i * stride] - wmean);
				W += wi ;
				wabsdev += (delta - wabsdev) * (wi / W);
			}
		}
		
		return wabsdev;		
	}
	
	double CStatis::wabsdev() const
	{
		const double wmean = get_wmean();
		return wabsdev(wmean);
	}
	
	double CStatis::wskew( const double& wmean, const double wsd ) const
	{
		/* Compute the weighted skewness of a dataset */
		
		double wskew = 0;
		double W = 0;
		
		size_t i;
		
		/* find the sum of the cubed deviations, normalized by the sd. */
		
		/* we use a recurrence relation to stably update a running value so
		there aren't any large sums that can overflow */
		
		for (i = 0; i < size; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0) {
				const double x = (data[i * stride] - wmean) / wsd;
				W += wi ;
				wskew += (x * x * x - wskew) * (wi / W);
			}
		}
		
		return wskew;		
	}
	
	double CStatis::wskew() const
	{
		const double wmean = get_wmean();
		const double wsd_ = wsd(wmean);
		return wskew(wmean, wsd_);
	}
	
	double CStatis::wkurtosis( const double& wmean, const double wsd ) const
	{
		/* takes a dataset and finds the kurtosis */
		
		long double wavg = 0, kurtosis;
		long double W = 0;
		size_t i;
		
		/* find the fourth moment the deviations, normalized by the sd */
		
		/* we use a recurrence relation to stably update a running value so
		there aren't any large sums that can overflow */
		
		for (i = 0; i < size; i++)
		{
			double wi = wdata[i * wstride];
			
			if (wi > 0) {
				const long double x = (data[i * stride] - wmean) / wsd;
				W += wi ;
				wavg += (x * x * x * x - wavg) * (wi / W);
			}
		}
		
		kurtosis = wavg - 3.0;  /* makes kurtosis zero for a Gaussian */
		
		return kurtosis;		
	}
	
	double CStatis::wkurtosis() const
	{
		const double wmean = get_wmean();
		const double wsd_ = wsd(wmean);
		return wkurtosis(wmean, wsd_);
	}

	
	double CStatis::max(size_t* mIndex /* = NULL*/) const
	{
		/* finds the index of the largest member of a dataset */
		/* if there is more than one largest value then we choose the first */
		
		double max = data[0 * stride];
		size_t i, max_index = 0;
		
		for (i = 0; i < size; i++)
		{
			double xi = data[i * stride];
			
			if (xi > max)
			{
				max = xi;
				max_index = i;
			}
			
		#ifdef FP
			if (isnan (xi))
			{
				return i;
			}
		#endif
		}
		if (mIndex)
		{
			*mIndex = max_index;
		}
		return max;
	}
	
	double CStatis::min(size_t* mIndex /*= NULL*/) const
	{
		/* finds the index of the smallest member of a dataset */
		/* if there is more than one largest value then we choose the first  */
		
		double min = data[0 * stride];
		size_t i, min_index = 0;
		
		for (i = 0; i < size; i++)
		{
			double xi = data[i * stride];
			
			if (xi < min)
			{
				min = xi;
				min_index = i;
			}
			
		#ifdef FP
			if (isnan (xi))
			{
				return i;
			}
		#endif
		}
		if (mIndex)
		{
			*mIndex = min_index;
		}
		return min;	
	}
	
	void CStatis::minmax( double & min_out, double& max_out, size_t* min_index_out /*= NULL*/, size_t* max_index_out/*= NULL*/ ) const
	{
		/* finds the smallest and largest members of a dataset */
		
		double min = data[0 * stride];
		double max = data[0 * stride];
		size_t i, min_index = 0, max_index = 0;
		
		for (i = 0; i < size; i++)
		{
			double xi = data[i * stride];
			
			if (xi < min)
			{
				min = xi;
				min_index = i;
			}
			
			if (xi > max)
			{
				max = xi;
				max_index = i;
			}
			
		#ifdef FP
			if (isnan (xi))
			{
				min_index = i;
				max_index = i;
				break;
			}
		#endif
		}
		if (min_index_out && max_index_out)
		{
			*min_index_out = min_index;
			*max_index_out = max_index;
		}
		min_out = min;
		max_out = max;
	}
}




