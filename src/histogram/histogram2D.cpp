/********************************************************************
	filename: 	histogram2D.cpp
	author:		huzhijian
	created:	5:5:2010   16:44
	brief:	2Î¬Öù×´Í¼
*********************************************************************/

#include <math.h>
#include <gsl_errno.h> 
#include <histogram2D.h>

namespace gslcpp
{
	CHistogram2D::CHistogram2D(const size_t nx, const size_t ny,bool b_zero)
		:gsldata(NULL), gslpdf(NULL)
	{
		gsldata = (b_zero == false) ? Alloc(nx, ny) : Calloc(nx, ny);
	}
	
	CHistogram2D::CHistogram2D( const size_t nx, const size_t ny,
								const double xmin, const double xmax,
                                const double ymin, const double ymax)
		:gsldata(NULL), gslpdf(NULL)
	{
		
		
		if (xmin >= xmax)
		{
			gsl_error("xmin must be less than xmax", __FILE__, __LINE__, GSL_EINVAL);
		}
		
		if (ymin >= ymax)
		{
			gsl_error ("ymin must be less than ymax", __FILE__, __LINE__, GSL_EINVAL);
		}
		
		gsldata = Calloc(nx, ny);
		
		if (gsldata == 0)
		{
			gsl_error("memory allocation failed", __FILE__, __LINE__, GSL_FAILURE);
		}
		
		make_uniform (gsldata->xrange, nx, xmin, xmax);
		make_uniform (gsldata->yrange, ny, ymin, ymax);
		
		return;
	}
	
	CHistogram2D::gsl_histogram2d* CHistogram2D::Alloc(const size_t nx, const size_t ny)
	{
		gsl_histogram2d *h;
		
		if (nx == 0)
		{
			GSL_ERROR_VAL ("histogram2d length nx must be positive integer",
				GSL_EDOM, 0);
		}
		
		if (ny == 0)
		{
			GSL_ERROR_VAL ("histogram2d length ny must be positive integer",
				GSL_EDOM, 0);
		}
		
		h = (gsl_histogram2d *) new gsl_histogram2d;
		
		if (h == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for histogram2d struct",
				GSL_ENOMEM, 0);
		}
		
		h->xrange = (double *) new double[nx+1];
		
		if (h->xrange == 0)
		{
			delete h;         /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d x ranges",
				GSL_ENOMEM, 0);
		}
		
		h->yrange = (double *) new double[ny+1];
		
		if (h->yrange == 0)
		{
			delete [] h->xrange;
			delete h;         /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d y ranges",
				GSL_ENOMEM, 0);
		}
		
		h->bin = new double[nx*ny];
		
		if (h->bin == 0)
		{
			delete [] h->xrange;
			delete [] h->yrange;
			delete h;        /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram bins",
				GSL_ENOMEM, 0);
		}
		
		h->nx = nx;
		h->ny = ny;
		
		return h;
	}
	
	
	void CHistogram2D::Free()
	{
		if (gsldata)
		{
			if (gsldata->xrange)
			{
				delete [] gsldata->xrange;
			}
			if (gsldata->yrange)
			{
				delete [] gsldata->yrange;
			}
			if (gsldata->bin)
			{
				delete[] gsldata->bin;
			}
			delete gsldata;
			gsldata = NULL;
		}

		if (gslpdf)
		{
			if (gslpdf->xrange)
			{
				delete [] gslpdf->xrange;
			}
			if (gslpdf->yrange)
			{
				delete [] gslpdf->yrange;
			}
			if (gslpdf->sum)
			{
				delete [] gslpdf->sum;
			}
			delete gslpdf;
			gslpdf = NULL;
		}
			
	}
	
	int CHistogram2D::SetRanges( const double xrange[], size_t xsize,
								 const double yrange[], size_t ysize)
	{
		size_t i;
		const size_t nx = gsldata->nx, ny = gsldata->ny;
		
		if (xsize != (nx + 1))
		{
			GSL_ERROR_VAL ("size of xrange must match size of histogram", 
				GSL_EINVAL, 0);
		}
		
		if (ysize != (ny + 1))
		{
			GSL_ERROR_VAL ("size of yrange must match size of histogram", 
				GSL_EINVAL, 0);
		}
		
		/* initialize ranges */
		
		for (i = 0; i <= nx; i++)
		{
			gsldata->xrange[i] = xrange[i];
		}
		
		for (i = 0; i <= ny; i++)
		{
			gsldata->yrange[i] = yrange[i];
		}
		
		/* clear contents */
		
		for (i = 0; i < nx * ny; i++)
		{
			gsldata->bin[i] = 0;
		}
		
		return GSL_SUCCESS;
	}
	
	int CHistogram2D::SetRanges( double xmin, double xmax,
                                 double ymin, double ymax)
	{
		size_t i;
		const size_t nx = gsldata->nx, ny = gsldata->ny;
		
		if (xmin >= xmax)
		{
			GSL_ERROR_VAL ("xmin must be less than xmax", GSL_EINVAL, 0);
		}
		
		if (ymin >= ymax)
		{
			GSL_ERROR_VAL ("ymin must be less than ymax", GSL_EINVAL, 0);
		}
		
		/* initialize ranges */
		
		make_uniform (gsldata->xrange, nx, xmin, xmax);
		make_uniform (gsldata->yrange, ny, ymin, ymax);
		
		/* clear contents */
		
		for (i = 0; i < nx * ny; i++)
		{
			gsldata->bin[i] = 0;
		}
		
		return GSL_SUCCESS;
	}
	
	void CHistogram2D::make_uniform(double range[], const size_t n,const double xmin, const double xmax)
	{
		size_t i;
		
		for (i = 0; i <= n; i++)
		{
			double f1 = ((double) (n-i) / (double) n);
			double f2 = ((double) i / (double) n);
			range[i] = f1 * xmin +  f2 * xmax;
		}		
	}
	
	CHistogram2D::gsl_histogram2d* CHistogram2D::Calloc( const size_t nx, const size_t ny )
	{
		gsl_histogram2d *h;
		
		if (nx == 0)
		{
			GSL_ERROR_VAL ("histogram2d length nx must be positive integer",
				GSL_EDOM, 0);
		}
		
		if (ny == 0)
		{
			GSL_ERROR_VAL ("histogram2d length ny must be positive integer",
				GSL_EDOM, 0);
		}
		
		h = (gsl_histogram2d *) new gsl_histogram2d;
		
		if (h == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for histogram2d struct",
				GSL_ENOMEM, 0);
		}
		
		h->xrange = (double *) new double[nx+1];
		
		if (h->xrange == 0)
		{
			delete h;         /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d x ranges",
				GSL_ENOMEM, 0);
		}
		
		h->yrange = (double *) new double[ny+1];
		
		if (h->yrange == 0)
		{
			delete [] h->xrange;
			delete h;         /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d y ranges",
				GSL_ENOMEM, 0);
		}
		
		h->bin = (double *) new double[nx*ny];
		
		if (h->bin == 0)
		{
			delete [] h->yrange;
			delete [] h->xrange;
			delete h;         /* exception in constructor, avoid memory leak */
			h = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram bins",
				GSL_ENOMEM, 0);
		}
		
		{
			size_t i;
			
			for (i = 0; i < nx + 1; i++)
			{
				h->xrange[i] = i;
			}
			
			for (i = 0; i < ny + 1; i++)
			{
				h->yrange[i] = i;
			}
			
			for (i = 0; i < nx * ny; i++)
			{
				h->bin[i] = 0;
			}
		}
		
		h->nx = nx;
		h->ny = ny;
		
		return h;
	}
	
	void CHistogram2D::Resize( const size_t nx, const size_t ny )
	{
		if (nx == gsldata->nx && ny == gsldata->ny)
		{
			return;
		}
		Free();
		gsldata = Alloc(nx, ny);
	}
	
	CHistogram2D& CHistogram2D::operator=( const CHistogram2D & src )
	{
		const size_t nx = src.gsldata->nx;
		const size_t ny = src.gsldata->ny;
		const size_t n = nx*ny;
		Resize(nx, ny);
		size_t i;
		for (i = 0; i <= nx; i++)
		{
			gsldata->xrange[i] = src.gsldata->xrange[i];
		}
		
		for (i = 0; i <= ny; i++)
		{
			gsldata->yrange[i] = src.gsldata->yrange[i];
		}

		for (i = 0; i < n; i++)
		{
			gsldata->bin[i] = src.gsldata->bin[i];
		}
		
		return *this;		
	}
	
	double CHistogram2D::Get( const size_t i, const size_t j) const
	{
		
		if (i >= gsldata->nx)
		{
			GSL_ERROR_VAL ("index i lies outside valid range of 0 .. nx - 1",
				GSL_EDOM, 0);
		}
		
		if (j >= gsldata->ny)
		{
			GSL_ERROR_VAL ("index j lies outside valid range of 0 .. ny - 1",
				GSL_EDOM, 0);
		}
		
		return gsldata->bin[i * gsldata->ny + j];
	}
	
	int CHistogram2D::GetXRange( const size_t index, double &xlower, double&xupper ) const
	{
		if (index >= gsldata->nx)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. nx - 1", GSL_EDOM);
		}
		
		xlower = gsldata->xrange[index];
		xupper = gsldata->xrange[index + 1];
		
		return GSL_SUCCESS;		
	}
	
	int CHistogram2D::GetYRange( const size_t index, double &ylower, double&yupper ) const
	{
		if (index >= gsldata->ny)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. ny - 1", GSL_EDOM);
		}
		
		ylower = gsldata->yrange[index];
		yupper = gsldata->yrange[index + 1];
		
		return GSL_SUCCESS;		
	}

	size_t CHistogram2D::Find(const size_t n, const double range[], const double x ) const
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
	inline double CHistogram2D::xMax() const
	{
		return gsldata->xrange[gsldata->nx];	
	}

	inline double CHistogram2D::xMin() const
	{
		return gsldata->xrange[0];
	}
	
	inline double CHistogram2D::yMax() const
	{
		return gsldata->yrange[gsldata->ny];	
	}
	
	inline double CHistogram2D::yMin() const
	{
		return gsldata->yrange[0];
	}

	inline double CHistogram2D::nX() const
	{
		return gsldata->nx;		
	}
	
	inline double CHistogram2D::nY() const
	{
		return gsldata->ny;		
	}

	void CHistogram2D::Reset()
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		double* bin = gsldata->bin;
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] = 0;
		}
	}
	
	int CHistogram2D::Increment( const double x, const double y )
	{
		return Increment (x, y, 1.0);	
	}
	
	int CHistogram2D::Increment( const double x, const double y, const double weight)
	{
		size_t i = 0, j = 0;
		i = Find(gsldata->nx, gsldata->xrange, x);
		j = Find(gsldata->ny, gsldata->yrange, y);
		
		if (i >= gsldata->nx)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. nx - 1",
				GSL_ESANITY);
		}
		
		if (j >= gsldata->ny)
		{
			GSL_ERROR ("index lies outside valid range of 0 .. ny - 1",
				GSL_ESANITY);
		}

		gsldata->bin[i * gsldata->ny + j] += weight;

		return GSL_SUCCESS;		
	}
	
	double CHistogram2D::MaxBin( size_t *iIndex /*= NULL*/ , size_t* jIndex /*= NULL*/) const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t imax = 0, jmax = 0, i, j;
		double max = gsldata->bin[0 * ny + 0];
		
		for (i = 0; i < nx; i++)
		{
			for (j = 0; j < ny; j++)
			{
				double x = gsldata->bin[i * ny + j];
				
				if (x > max)
				{
					max = x;
					imax = i;
					jmax = j;
				}
			}
		}
		if (iIndex != NULL && jIndex != NULL)
		{
			*iIndex = imax;
			*jIndex = jmax;	
		}
		
		return max;
	}
	
	double CHistogram2D::MinBin( size_t *iIndex /*= NULL*/ , size_t* jIndex /*= NULL*/) const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t imin = 0, jmin = 0, i, j;
		double min = gsldata->bin[0 * ny + 0];
		
		for (i = 0; i < nx; i++)
		{
			for (j = 0; j < ny; j++)
			{
				double x = gsldata->bin[i * ny + j];
				
				if (x < min)
				{
					min = x;
					imin = i;
					jmin = j;
				}
			}
		}
		if (iIndex != NULL && jIndex != NULL)
		{
			*iIndex = imin;
			*jIndex = jmin;	
		}
		return min;	
	}
	
	double CHistogram2D::Sum() const
	{
		const size_t n = gsldata->nx * gsldata->ny;
		double sum = 0;
		size_t i = 0;
		
		while (i < n)
			sum += gsldata->bin[i++];
		
		return sum;	
	}
	
	double CHistogram2D::xMean() const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t i;
		size_t j;
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wmean = 0;
		long double W = 0;
		const double* xrange = gsldata->xrange;
		const double* bin = gsldata->bin;
		for (i = 0; i < nx; i++)
		{
			double xi = (xrange[i + 1] + xrange[i]) / 2.0;
			double wi = 0;
			
			for (j = 0; j < ny; j++)
			{
				double wij = bin[i * ny + j];
				if (wij > 0)
					wi += wij;
			}
			if (wi > 0)
			{
				W += wi;
				wmean += (xi - wmean) * (wi / W);
			}
		}
		
		return wmean;		
	}
	
	double CHistogram2D::yMean() const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t i;
		size_t j;
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wmean = 0;
		long double W = 0;
		const double* yrange = gsldata->yrange;
		const double* bin = gsldata->bin;
		for (j = 0; j < ny; j++)
		{
			double yj = (yrange[j + 1] + yrange[j]) / 2.0;
			double wj = 0;
			
			for (i = 0; i < nx; i++)
			{
				double wij = bin[i * ny + j];
				if (wij > 0)
					wj += wij;
			}
			
			if (wj > 0)
			{
				W += wj;
				wmean += (yj - wmean) * (wj / W);
			}
		}
		
		return wmean;
	}

	double CHistogram2D::xSigma() const
	{
		const double xmean = xMean();
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t i;
		size_t j;
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wvariance = 0;
		long double W = 0;
		const double* xrange = gsldata->xrange;
		const double* bin = gsldata->bin;

		for (i = 0; i < nx; i++)
		{
			double xi = (xrange[i + 1] + xrange[i]) / 2 - xmean;
			double wi = 0;
			
			for (j = 0; j < ny; j++)
			{
				double wij = bin[i * ny + j];
				if (wij > 0)
					wi += wij;
			}
			
			if (wi > 0)
			{
				W += wi;
				wvariance += ((xi * xi) - wvariance) * (wi / W);
			}
		}
		
		{
			double xsigma = sqrt (wvariance);
			return xsigma;
		}
	}
	
	double CHistogram2D::ySigma() const
	{
		const double ymean = yMean ();
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t i;
		size_t j;
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wvariance = 0;
		long double W = 0;
		const double* yrange = gsldata->yrange;
		const double* bin = gsldata->bin;

		for (j = 0; j < ny; j++)
		{
			double yj = (yrange[j + 1] + yrange[j]) / 2.0 - ymean;
			double wj = 0;
			
			for (i = 0; i < nx; i++)
			{
				double wij = bin[i * ny + j];
				if (wij > 0)
					wj += wij;
			}
			if (wj > 0)
			{
				W += wj;
				wvariance += ((yj * yj) - wvariance) * (wj / W);
			}
		}
		
		{
			double ysigma = sqrt (wvariance);
			return ysigma;
		}
	}
	bool CHistogram2D::operator==( const CHistogram2D& other ) const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		if ((nx != other.gsldata->nx) || (ny != other.gsldata->ny))
		{
			return false;
		}
		{
			size_t i;
			const double* xrange = gsldata->xrange;
			const double* yrange = gsldata->yrange;

			/* init ranges */
			for (i = 0; i <= (nx); i++)
			{
				if (xrange[i] != xrange[i])
				{
					return false;
				}
			}
			for (i = 0; i <= (ny); i++)
			{
				if (yrange[i] != yrange[i])
				{
					return false;
				}
			}
		}
		return true;		
	}
	
	bool CHistogram2D::operator!=( const CHistogram2D& other ) const
	{
		return !(*this == other);		
	}
	
	CHistogram2D& CHistogram2D::operator+=( const CHistogram2D & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;

		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] += bin_other[i];
		}
		
		return *this;
	}
	
	CHistogram2D& CHistogram2D::operator+=( const double & other )
	{
		double *bin = gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;
		
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] += other;
		}
		
		return *this;	
	}
	
	CHistogram2D& CHistogram2D::operator-=( const CHistogram2D & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;
		
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] -= bin_other[i];
		}
		
		return *this;	
	}
	
	CHistogram2D& CHistogram2D::operator*=( const CHistogram2D & other )
	{
		
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;
		
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] *= bin_other[i];
		}
		
		return *this;	
	}
	
	CHistogram2D& CHistogram2D::operator*=( const double & other )
	{
		double *bin = gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;
		
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] *= other;
		}
		
		return *this;
	}
	
	CHistogram2D& CHistogram2D::operator/=( const CHistogram2D & other )
	{
		if (*this != other)
		{
			gsl_error("histograms have different binning", __FILE__, __LINE__, GSL_EINVAL);
		}
		double *bin = gsldata->bin;
		const double *bin_other = other.gsldata->bin;
		const size_t nx  = gsldata->nx;
		const size_t ny  = gsldata->ny;
		
		for (size_t i = 0; i < nx*ny; i++)
		{
			bin[i] /= bin_other[i];
		}
		
		return *this;		
	}
	
	int CHistogram2D::Save( const char* filename ) const
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		const double* bin = gsldata->bin;
		const double* xrange = gsldata->xrange;
		const double* yrange = gsldata->yrange;
		FILE* fp = fopen(filename, "w");

		for (size_t i =0; i<nx; i++)
		{
			for (size_t j=0; j<ny; j++)
			{
				if (fprintf(fp, "%g\t%g\t%g\t%g\t%g\n", xrange[i], xrange[i+1], yrange[j], yrange[j+1],bin[i * ny + j]) < 0)
				{
					fclose(fp);
					GSL_ERROR ("Histogram save failed", GSL_EFAILED);
				}
			}
			
		}
		fclose(fp);
		return GSL_SUCCESS;		
	}
	
	int CHistogram2D::Load( const char* filename )
	{
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		double* bin = gsldata->bin;
		double* xrange = gsldata->xrange;
		double* yrange = gsldata->yrange;
		double xupper, yupper;
		FILE* fp = fopen(filename, "r");
		
		for (size_t i =0; i<nx; i++)
		{
			for (size_t j=0; j<ny; j++)
			{
				if (fscanf(fp, "%lg %lg %lg %lg %lg\n", 
					xrange + i, &xupper,
					yrange + j, &yupper, 
					bin + i * ny + j) != 5)
				{
					fclose(fp);
					GSL_ERROR ("Histogram load failed", GSL_EFAILED);
				}
			}
			yrange[ny] = yupper;
		}
			
		
		fclose(fp);
		xrange[nx] = xupper;

		return GSL_SUCCESS;
	}
	
	CHistogram2D::gsl_histogram2d_pdf* CHistogram2D::AllocPdf(const size_t nx, const size_t ny)
	{
		const size_t n = nx * ny;
		
		gsl_histogram2d_pdf *p;
		
		if (n == 0)
		{
			GSL_ERROR_VAL ("histogram2d pdf length n must be positive integer",
				GSL_EDOM, 0);
		}
		
		p = (gsl_histogram2d_pdf *) new gsl_histogram2d_pdf;
		
		if (p == 0)
		{
			GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf struct",
				GSL_ENOMEM, 0);
		}
		
		p->xrange = (double *) new double[nx+1];
		
		if (p->xrange == 0)
		{
			delete p;         /* exception in constructor, avoid memory leak */
			p = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf xranges",
				GSL_ENOMEM, 0);
		}
		
		p->yrange = (double *) new double[ny+1];
		
		if (p->yrange == 0)
		{
			delete [] p->xrange;
			delete p;        /* exception in constructor, avoid memory leak */
			p = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf yranges",
				GSL_ENOMEM, 0);
		}
		
		p->sum = (double *) new double[n+1];
		
		if (p->sum == 0)
		{
			delete [] p->xrange;
			delete [] p->yrange;
			delete p;        /* exception in constructor, avoid memory leak */
			p = NULL;
			GSL_ERROR_VAL ("failed to allocate space for histogram2d pdf sums",
				GSL_ENOMEM, 0);
		}
		
		p->nx = nx;
		p->ny = ny;
		
		return p;
	}
	
	int CHistogram2D::InitPdf()
	{
		size_t i;
		const double* bin = gsldata->bin;
		const double* xrange = gsldata->xrange;
		const double* yrange = gsldata->yrange;
		double* xrange_pdf = gslpdf->xrange;
		double* yrange_pdf = gslpdf->yrange;
		double* sum_pdf = gslpdf->sum;
		const size_t nx = gslpdf->nx;
		const size_t ny = gslpdf->ny;
		const size_t n = nx * ny;
		
		if (nx != gsldata->nx || ny != gsldata->ny)
		{
			GSL_ERROR ("histogram2d size must match pdf size", GSL_EDOM);
		}
		
		for (i = 0; i < n; i++)
		{
			if (bin[i] < 0)
			{
				GSL_ERROR ("histogram bins must be non-negative to compute"
					"a probability distribution", GSL_EDOM);
			}
		}
		
		for (i = 0; i < nx + 1; i++)
		{
			xrange_pdf[i] = xrange[i];
		}
		
		for (i = 0; i < ny + 1; i++)
		{
			yrange_pdf[i] = yrange[i];
		}
		
		{
			double mean = 0, sum = 0;
			
			for (i = 0; i < n; i++)
			{
				mean += (bin[i] - mean) / ((double) (i + 1));
			}
			
			sum_pdf[0] = 0;
			
			for (i = 0; i < n; i++)
			{
				sum += (bin[i] / mean) / n;
				sum_pdf[i + 1] = sum;
			}
		}
		
		return GSL_SUCCESS;	
	}
	
	int CHistogram2D::PdfSample( double r1, double r2, double &x, double &y)
	{
		if (!gslpdf)
		{
			gslpdf = AllocPdf(gsldata->nx, gsldata->ny);
			InitPdf();
		}
		size_t k;
		
		/* Wrap the exclusive top of the bin down to the inclusive bottom of
		the bin. Since this is a single point it should not affect the
		distribution. */
		
		if (r2 == 1.0)
		{
			r2 = 0.0;
		}
		if (r1 == 1.0)
		{
			r1 = 0.0;
		}
		
		k = Find(gslpdf->nx * gslpdf->ny, gslpdf->sum, r1);
		
		
		size_t i = k / gslpdf->ny;
		size_t j = k - (i * gslpdf->ny);
		double delta = (r1 - gslpdf->sum[k]) / (gslpdf->sum[k + 1] - gslpdf->sum[k]);
		x = gslpdf->xrange[i] + delta * (gslpdf->xrange[i + 1] - gslpdf->xrange[i]);
		y = gslpdf->yrange[j] + r2 * (gslpdf->yrange[j + 1] - gslpdf->yrange[j]);
		return GSL_SUCCESS;		
	}
	
	double CHistogram2D::Cov() const
	{
		const double xmean = xMean ();
		const double ymean = yMean ();
		const size_t nx = gsldata->nx;
		const size_t ny = gsldata->ny;
		size_t i;
		size_t j;
		
		/* Compute the bin-weighted arithmetic mean M of a histogram using the
		recurrence relation
		
		  M(n) = M(n-1) + (x[n] - M(n-1)) (w(n)/(W(n-1) + w(n))) 
		  W(n) = W(n-1) + w(n)
		  
		*/
		
		long double wcovariance = 0;
		long double W = 0;
		const double* xrange = gsldata->xrange;
		const double* yrange = gsldata->yrange;
		const double* bin = gsldata->bin;
		for (j = 0; j < ny; j++)
		{
			for (i = 0; i < nx; i++)
			{
				double xi = (xrange[i + 1] + xrange[i]) / 2.0 - xmean;
				double yj = (yrange[j + 1] + yrange[j]) / 2.0 - ymean;
				double wij = bin[i * ny + j];
				
				if (wij > 0)
				{
					W += wij;
					wcovariance += ((xi * yj) - wcovariance) * (wij / W);
				}
			}
		}
		
		return wcovariance;		
	}
	
	CHistogram2D::~CHistogram2D()
	{
		Free();		
	}
}