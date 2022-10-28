/********************************************************************
	filename: 	Statis.h
	author:		hu zhijian
	created:	5:5:2010   16:43
	brief:	统计类
*********************************************************************/

#ifndef NSL_STATIS_H__
#define NSL_STATIS_H__


#include <NSL.h>

#include <stdlib.h>

namespace gslcpp
{
	class NSL_EXPORT CStatis//: public CVectorView
	{
	private:
		const double* data;
		const size_t stride;
		const size_t size;

		double* wdata;
		size_t wstride;
		size_t wsize;
		double compute_variance(double mean) const;
		double compute_covariance(const CStatis& other , const double& mean1, const double& mean2 )const;
		double compute_wvariance(const double& wmean) const;
		double compute_factor() const;
	public:
		CStatis(const double * base, size_t n, size_t stride = 1);
		virtual ~CStatis();

		void AddWight(const double * base, size_t n, size_t stride = 1);
		/************************************************************************/
		/* 平均值                                                               
		/************************************************************************/
		double get_mean() const;
		
		/************************************************************************/
		/* 方差                                                                    
		/************************************************************************/
		double variance() const ;
		double variance(double mean) const ;
		double douvariance_with_fixed_mean(const double & mean) const;
		/************************************************************************/
		/* standard deviation                                                        
		/************************************************************************/
		double std_deviation() const;
		double std_deviation(const double& mean) const;
		double std_deviation_with_fixed_mean(const double& mean) const;

		/************************************************************************/
		/* absolute deviation                                                            
		/************************************************************************/
		double absdev() const;
		double absdev(const double& mean) const;

		/************************************************************************/
		/* Higher moments (skewness and kurtosis)                                               
		/************************************************************************/
		double skew() const;
		double skew(const double& mean, const double& sd) const;
		double kurtosis() const;
		double kurtosis(const double& mean, const double& sd) const;

		/************************************************************************/
		/* Autocorrelation                                       
		/************************************************************************/
		double lag1_autocorrelation() const;
		double lag1_autocorrelation(const double& mean) const;

		/************************************************************************/
		/* Covariance                                                                
		/************************************************************************/
		double covariance(const CStatis& other) const;
		double covariance(const CStatis& other, const double& mean1, const double& mean2) const;

		/************************************************************************/
		/* Weighted Samples                                                           
		/************************************************************************/
		double get_wmean() const;
		double wvariance() const;
		double wvariance(const double& wmean) const;
		double wsd() const;
		double wsd(const double& wmean) const;
		double wvariance_with_fixed_mean(const double& wmean) const;
		double wsd_with_fixed_mean(const double& wmean) const;
		double wabsdev() const;
		double wabsdev(const double& wmean) const;
		double wskew() const;
		double wskew(const double& wmean, const double wsd) const;
		double wkurtosis() const;
		double wkurtosis(const double& wmean, const double wsd) const;

		/************************************************************************/
		/* Maximum and Minimum values                                       
		/************************************************************************/
		double max(size_t* mIndex = NULL) const;
		double min(size_t* mIndex = NULL) const;
		void minmax(double & min_out, double& max_out, size_t* min_index_out = NULL, size_t* max_index_out= NULL) const;
	};


}


#endif // NSL_STATIS_H__