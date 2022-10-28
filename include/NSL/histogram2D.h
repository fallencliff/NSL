/********************************************************************
	filename: 	histogram2D.h
	author:		hu zhijian
	created:	5:5:2010   16:44
	brief:	2维柱状图
*********************************************************************/

#ifndef NSL_HISTOGRAM2D_H__
#define NSL_HISTOGRAM2D_H__

#include <NSL.h>

namespace gslcpp
{
	class NSL_EXPORT CHistogram2D
	{
	public:
		typedef struct 
		{
			size_t nx, ny ;
			double * xrange ;
			double * yrange ;
			double * bin ;
		} gsl_histogram2d;
		
		typedef struct 
		{
			size_t nx, ny ;
			double * xrange ;
			double * yrange ;
			double * sum ;
		} gsl_histogram2d_pdf;

	private:
		gsl_histogram2d* gsldata;
		gsl_histogram2d_pdf* gslpdf;

		gsl_histogram2d* Alloc(const size_t nx, const size_t ny);
		gsl_histogram2d* Calloc(const size_t nx, const size_t ny);
		gsl_histogram2d_pdf* AllocPdf(const size_t nx, const size_t ny);
		int InitPdf();
		void Free();
		void make_uniform(double range[], const size_t n,const double xmin, const double xmax);
		void Resize(const size_t nx, const size_t ny );
	public:
		CHistogram2D(const size_t nx , const size_t ny, bool b_zero = false);
		CHistogram2D( const size_t nx, const size_t ny,
					  const double xmin, const double xmax,
                      const double ymin, const double ymax);

		virtual ~CHistogram2D();

		int SetRanges(const double xrange[], size_t xsize,
					  const double yrange[], size_t ysize);
		int SetRanges(double xmin, double xmax,
                      double ymin, double ymax);

		size_t Find(const size_t n, const double range[], const double x) const;

		double Get(const size_t i, const size_t j) const;
		int GetXRange(const size_t index, double &xlower, double&xupper)const;
		int GetYRange(const size_t index, double &ylower, double&yupper) const;
		double xMax() const;
		double xMin() const;
		double yMax() const;
		double yMin() const;
		double nX() const;
		double nY() const;
		
		void Reset();

		int Increment(const double x, const double y);
		int Increment(const double x, const double y, const double weight);

		double MaxBin(size_t *iIndex = NULL, size_t* jIndex = NULL) const;
		double MinBin(size_t *iIndex = NULL, size_t* jIndex = NULL) const;

		/************************************************************************/
		/* 操作符重载
		/************************************************************************/
		CHistogram2D& operator =(const CHistogram2D & src);
		bool operator ==(const CHistogram2D& other) const; //比较n和range
		bool operator !=(const CHistogram2D& other) const;
		CHistogram2D& operator +=(const CHistogram2D & other);
		CHistogram2D& operator +=(const double & other);
		CHistogram2D& operator -=(const CHistogram2D & other);
		CHistogram2D& operator *=(const CHistogram2D & other);
		CHistogram2D& operator *=(const double & other);
		CHistogram2D& operator /=(const CHistogram2D & other);

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
		double xMean() const;
		double yMean() const;

		/************************************************************************/
		/*This function returns the standard deviation of the histogrammed variable, where the
		/*histogram is regarded as a probability distribution                                                                      
		/************************************************************************/
		double xSigma() const;
		double ySigma() const;

		int PdfSample(double r1, double r2, double &x, double &y);

		double Cov() const ;
	};
}



#endif // NSL_HISTOGRAM2D_H__