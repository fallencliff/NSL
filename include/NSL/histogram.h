/********************************************************************
	filename: 	histogram.h
	author:		hu zhijian
	created:	5:5:2010   16:43
	brief:	柱状图
*********************************************************************/

#ifndef NSL_HISTOGRAM_H__
#define NSL_HISTOGRAM_H__


#include <SL_dll.h>
#include <CommonStruct.h>

namespace gslcpp
{
	class NSL_EXPORT CHistogram
	{
	public:

	private:
		gsl_histogram* gsldata;
		gsl_histogram_pdf* gslpdf;

		gsl_histogram* Alloc(const size_t n);
		gsl_histogram* Calloc(const size_t n);
		gsl_histogram_pdf* AllocPdf(const size_t n);
		int InitPdf();
		void Free();
		void make_uniform(const size_t n,const double xmin, const double xmax);
		void Resize(const size_t n);
	public:
		CHistogram(const size_t n , bool b_zero = false);
		CHistogram(const size_t n, const double xmin, const double xmax);
		virtual ~CHistogram();

		int SetRanges(const double range[], const size_t size);
		int SetRanges(const double xmin, const double xmax);

		size_t Find(const size_t n, const double range[], const double x) const;

		double Get(const size_t index) const;
		int GetRange(const size_t index, double &lower, double&upper )const;

		double Max() const;
		double Min() const;
		double Bins() const;
		
		void Reset();

		int Increment(const double x);
		int Increment(const double x, const double weight);

		/************************************************************************/
		/* gsl_histogram_max_val ,gsl_histogram_max_bin  
		/************************************************************************/
		double MaxBin(size_t *index = NULL) const;
		double MinBin(size_t *index = NULL) const;

		/************************************************************************/
		/* 操作符重载
		/************************************************************************/
		CHistogram& operator =(const CHistogram & src);
		bool operator ==(const CHistogram& other) const; //比较n和range
		bool operator !=(const CHistogram& other) const;
		CHistogram& operator +=(const CHistogram & other);
		CHistogram& operator +=(const double & other);
		CHistogram& operator -=(const CHistogram & other);
		CHistogram& operator *=(const CHistogram & other);
		CHistogram& operator *=(const double & other);
		CHistogram& operator /=(const CHistogram & other);

		int Save(const char* filename) const;
		int Load(const char* filename);

		/************************************************************************/
		/*  sum up all bins of histogram                                         */
		/************************************************************************/
		double Sum() const;

		/************************************************************************/
		/* This function returns the mean of the histogrammed variable, where the histogram
		/*	is regarded as a probability distribution 
		/************************************************************************/
		double Mean() const;

		/************************************************************************/
		/*This function returns the standard deviation of the histogrammed variable, where the
		/*histogram is regarded as a probability distribution                                                                      
		/************************************************************************/
		double Sigma() const;

		double PdfSample(double r);
	};
}


#endif // NSL_HISTOGRAM_H__