/********************************************************************
	filename: 	Matrix.h
	author:		hu zhijian
	created:	5:5:2010   16:45
	brief:	矩阵类声明
*********************************************************************/
#ifndef NSL_MATRIX_H__
#define NSL_MATRIX_H__

#include <time.h>
#include <gsl_errno.h>
#include <NSL.h>
#include <Vector.h>


namespace gslcpp
{
using gslcpp::VectorView;
using gslcpp::ConstVectorView;

class mslice 
	{
	public:
		mslice( const size_t startRow, const size_t startCol, const size_t rows, const size_t cols)
			: startRow_( startRow), startCol_( startCol),
			rows_( rows), cols_( cols) {}
		mslice( const mslice &ms)
			: startRow_( ms.startRow_), startCol_( ms.startCol_),
			rows_( ms.rows_), cols_( ms.cols_) {}
		
		size_t size() const {return rows_ * cols_;}
		size_t rows() const {return rows_;}
		size_t cols() const {return cols_;}
		size_t startRow() const {return startRow_;}
		size_t startCol() const {return startCol_;}
		size_t endRow() const {return startRow_ + rows_;}
		size_t endCol() const {return startCol_ + cols_;}
	private:
		size_t startRow_;
		size_t startCol_;
		size_t rows_;
		size_t cols_;
	};

typedef enum 
{
    CC_HORIZONTAL = 0,
	CC_VERTICAL   = 1
}Cc_Dir;

template<typename T>
struct gDecompositionInfo 
{
	gDecompositionInfo() : rank( 0), determinant( T(0)) {}
	gDecompositionInfo( size_t r, T d)  : rank( r), determinant( T(d)) {}
	size_t rank;
	T      determinant;
};

	//前向声明
	template<typename T>
		class ConstMatrixView;
	template<typename T>
		class AlmostMatrix;
	template<typename T>
		class MatrixView;
	template<typename T>
		class Matrix;
	//
	// 矩阵赋值
	//
	template<typename T>
	struct CommaAssignmentForMatrix 
	{
		CommaAssignmentForMatrix( AlmostMatrix<T> &m, const T &t)
		{
#if NSL_RANGE_CHECK
			if ( m.size() == 0) 
			{
				throw gsl_error( "CommaAssignmentForMatrix: Matrix has not been sized.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			m_    = &m;
			index = 0;
			m_->ref(index/m.cols(), index%m.cols()) = t;
		}
		
		CommaAssignmentForMatrix& operator, ( const T &t)
		{
			++index;
#if NSL_RANGE_CHECK
			if ( index >= m_->size()) 
			{
				throw gsl_error( "CommaAssignmentForMatrix::operator,(t): Assignment outside matrix range.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			m_->ref(index/m_->cols(), index%m_->cols()) = t;
			return *this;
		}
		
		AlmostMatrix<T> *m_;
		size_t index;
		};

	//一般矩阵基类
	template<typename T>
		class GneMatrix
	{
	protected:
		
		//单线程运行时，可作为缓存使用
#ifndef NSL_SUPPORT_MULTITHREAD
		static Matrix<T> cacheMatrix0_;
		static Matrix<T> cacheMatrix1_;
		static Matrix<T> cacheMatrix2_;
#endif

	public:
		GneMatrix() {}
		virtual ~GneMatrix() {}

		//派生类要重写的虚函数
		virtual const mslice& getMatrixSlice() 	const = 0;
		virtual const valarray<T>& constVal() 	const = 0;
		virtual const mslice& getOrgSlice()		const = 0;
		
		//派生类公有的函数, all are const functions
		size_t size() 				const 
		{ 
			return getMatrixSlice().size() ;
		}
		size_t startRow()			const {	return getMatrixSlice().startRow() ;}
		size_t startCol()			const {	return getMatrixSlice().startCol() ;}

		size_t rows() 				const { return getMatrixSlice().rows() ;}
		size_t cols() 				const { return getMatrixSlice().cols() ;}
		
		size_t typesize()		const { return sizeof(T) ;}
		inline T cref(size_t i, size_t j) const { return constVal()[index(i, j)] ;}
		
		
		
		//在valarray中的绝对位置
		inline size_t index(size_t i, size_t j) const  
		{ return (getMatrixSlice().startRow()+i)*getOrgSlice().cols() + (j+startCol());}
		
		
		
		valarray<T> toVal()		const 
		{ 
			valarray<T> cache(size());
			size_t i=0;
			size_t cacheIndex = 0;
			for(i=0; i<rows(); i++)
			{
				for(size_t j=0; j<cols(); j++)
				{
					cache[cacheIndex++] = cref(i,j);
				}
				
			}
			return cache;
		} 
		
		T sum() 				const 
		{ 
			T cache = 0.0;
			size_t i=0,j=0;
			size_t ros_ = rows();
			size_t cols_ = cols();
			for (i=0; i<ros_; i++)
			{
				for (j=0; j<cols_; j++)
				{
					cache += cref(i,j);
				}
				
			}

			return cache;
		}
		T sumAbs() 				const 
		{ 
			T cache = 0.0;
			size_t i=0,j=0;
			size_t ros_ = rows();
			size_t cols_ = cols();
			for (i=0; i<ros_; i++)
			{
				for (j=0; j<cols_; j++)
				{
					cache += nMath::abs(cref(i,j));
				}
				
			}
			return cache;
		}
	
 		T min(size_t* i_index_out = NULL, size_t * j_index_out = NULL)				const 
		{ 
			T cache = cref(0,0);
			size_t i=0,j=0;
			size_t ros_ = rows();
			size_t cols_ = cols();
			size_t iIndex =0;
			size_t jIndex = 0;
			for (i=0; i<ros_; i++)
			{
				for (j=0; j<cols_; j++)
				{
					if (cref(i,j) < cache)
					{
						cache = cref(i,j);
						iIndex = i;
						jIndex = j;

					}
				}
				
			}
			if (i_index_out && j_index_out)
			{
				*i_index_out = iIndex;
				*j_index_out = jIndex;
			}
			return cache;
		}
 		T max(size_t* i_index_out = NULL, size_t * j_index_out = NULL)				const 
		{ 
			T cache = cref(0,0);
			size_t i=0,j=0;
			size_t ros_ = rows();
			size_t cols_ = cols();
			size_t iIndex =0;
			size_t jIndex = 0;
			for (i=0; i<ros_; i++)
			{
				for (j=0; j<cols_; j++)
				{
					if (cref(i,j) > cache)
					{
						cache = cref(i,j);
						iIndex = i;
						jIndex = j;
						
					}
				}
				
			}
			if (i_index_out && j_index_out)
			{
				*i_index_out = iIndex;
				*j_index_out = jIndex;
			}
			return cache;
		}
 		T minAbs(size_t* i_index_out = NULL, size_t j_index_out = NULL)				const 
		{ 
			T cache = cref(0,0);
			size_t i=0,j=0;
			size_t ros_ = rows();
			size_t cols_ = cols();
			size_t iIndex =0;
			size_t jIndex = 0;
			for (i=0; i<ros_; i++)
			{
				for (j=0; j<cols_; j++)
				{
					if (nMath::abs(cref(i,j)) < cache)
					{
						cache = cref(i,j);
						iIndex = i;
						jIndex = j;	
					}
				}
				
			}
			if (i_index_out && j_index_out)
			{
				*i_index_out = iIndex;
				*j_index_out = jIndex;
			}
			return cache;
		}
 		T maxAbs(size_t* i_index_out = NULL, size_t j_index_out = NULL)				const 
		{ 
		
		}
		T average() 	const {return sum() / size();}
		T norm() 	const 
		{ 
			T t(0);
			for ( size_t i=0; i<rows(); ++i) 
			{
				for (size_t j=0; j<cols(); ++j)
				{
					t += cref(i,j) * cref(i,j);
				}
				
			}
			return ::sqrt(t);
		}

		T trace( int i = 0) const 
		{
			valarray<T> d = diag(i); 
			return d.sum();
		}

		valarray<T> diag( size_t i = 0) 	const 
		{ 
#if	NSL_RANGE_CHECK
			if ( i >= cols() || rows() + i < 1) 
			{
				throw gsl_error( "Matrix<T>::diag(i): Index i out of range.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			size_t start;
			size_t length;
			if ( i >=0) 
			{
				start  = i;
				length = cols() - i;
				if ( length > rows()) 
				{
					length = rows();
				}
			}
			else {
				start  = -i * cols();
				length = rows() + i;
				if ( length > cols()) {
					length = cols();
				}
			}
			return constVal()[ slice( start, length, cols() + 1)];
		}
		
		
		//Unary operators
		Matrix<T> operator+ () const {return *this;}
		Matrix<T> operator- () const 
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize(rows(), cols());
			for (size_t i=0; i<rows(); i++)
			{
				for (size_t j=0; j<cols(); j++)
				{
					cacheMatrix0_(i,j) = -cref(i,j);
				}
			}

			return cacheMatrix0_;
		}
		
		
		
		bool isSquare() const 
		{
			return rows()==cols() ;
		}

 		bool isZero()	const 
		{ 
			size_t rows_ = rows();
			size_t cols_ = cols();
			for (size_t i=0; i<rows_; i++)
			{
				for (size_t j=0; j<cols_; j++)
				{
					if (nMath::abs(cref(i,j)) > nMath::ZERO)
					{
						return false;
					}
				}
			}

			return true;
		}
 		bool hasSubMatrix(const mslice & s)	const 
		{ 
			if (s.rows() > rows()
				||s.cols()> cols()
				||s.startRow() >= rows()
				||s.startCol() >= cols()
				||s.endRow() > rows()
				||s.endCol() > cols())
			{
				return false;
			}

			return true;
		
		}
 		bool isUnit() 	const 
		{ 
			if ( !isDiagonal()) 
			{
				return false;
			}
			for ( size_t i = 0; i < rows(); ++i) 
			{
				if ( nMath::abs( cref(i,i) - T(1)) > nMath::ZERO) 
				{
					return false;
				}
			}
			return true;
		}
 		bool isSingular() 	const
		{
			if ( !isSquare() )  
			{
				return false;
			}
			return ( nMath::abs( det()) < nMath::ZERO);
		}
 		bool isDiagonal() 	const 
		{ 
			if (!isSquare()) 
			{
				return false;
			}
			
			for ( size_t i = 0; i < rows(); ++i) 
			{
				for ( size_t j = 0; j < cols(); ++j) 
				{
					if ( i != j && nMath::abs(cref(i,j)) > nMath::ZERO) 
					{
						return false;
					}
				}
			}
			return true;
		
		}
 		bool isScalar() const 
		{
			if ( !isDiagonal()) 
			{
				return false;
			}
			T t = cref(0,0);
			for ( size_t i = 1; i < rows(); ++i)
			{
				if ( nMath::abs( cref(i,i) - t) > nMath::ZERO)
				{
					return false;
				}
			}
			return true;
		}

 		bool isSymmetric() 	const 
		{ 
			if (!isSquare()) 
			{
				return false;
			}
			for ( size_t i = 0; i < rows(); ++i)
			{
				for ( size_t j = 0; j < cols(); ++j) 
				{
					if ( nMath::abs( cref(i,j) - cref(j,i)) > nMath::ZERO) 
					{
						return false;
					}
				}
			}
			return true;
		}

		//反对称矩阵
 		bool isSkewSymmetric() const 
		{ 
			if ( !isSquare()) 
			{
				return false;
			}
			for ( size_t j = 1; j < cols(); ++j) 
			{
				for ( size_t i = 0; i < j; ++i) 
				{
					if ( nMath::abs( cref(i,j) + cref(j,i)) > nMath::ZERO)
					{
						return false;
					}
				}
			}
			return true;
		}
 		bool isUpperTriangular() const 
		{
			if ( !isSquare()) 
			{
				return false;
			}
			for ( size_t i = 1; i < rows(); ++i) 
			{
				for ( size_t j = 0; j < i - 1; ++j) 
				{
					if ( nMath::abs( cref(i,j)) > nMath::ZERO) 
					{
						return false;
					}
				}
			}
			return true;
		}
 		bool isLowerTriangular() const 
		{
			if ( !isSquare()) 
			{
				return false;
			}
			for ( size_t j = 1; j < cols(); ++j) 
			{
				for ( size_t i = 0; i < j - 1; ++i) 
				{
					if ( nMath::abs(cref(i,j)) > nMath::ZERO) 
					{
						return false;
					}
				}
			}
			return true;
		}

		Matrix<T> operator! () const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = *this;
			cacheMatrix0_.inverse();
			return cacheMatrix0_;
		}
		//转置
		Matrix<T> operator ~ () const
		{
			return transpose();
		}

		Matrix<T> transpose () const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<double> cacheMatrix0_;
#endif
			cacheMatrix0_.resize(rows(), cols());
			
			for ( size_t i = 0; i < rows(); ++i)
			{
				for ( size_t j = 0; j < cols(); ++j) 
				{
					cacheMatrix0_(j,i) = cref(i,j);
				}
			}
			return cacheMatrix0_;
		}
 		bool isRowOrthogonal() 	const 
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = (*this) * ~(*this);
		
			return cacheMatrix0_.isUnit();
		}
 		bool isColumnOrthogonal() const 
		{ 
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = ~(*this) * (*this);
			
			return cacheMatrix0_.isUnit();
		}

		//
		// Boolean operators.
		//
		friend bool operator== ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
			if ( (m1.rows() != m2.rows()) || (m1.cols() != m2.cols()))
			{
				return false;
			}
			for ( size_t i = 0; i < m1.rows(); ++i) 
			{
				for ( size_t j = 0; j < m1.cols(); ++j) 
				{
					if ( nMath::abs( m1.cref(i,j) - m2.cref(i,j)) > nMath::ZERO) 
					{
						return false;
					}
				}
			}
			return true;
		}
		
		friend bool operator!= ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
			return !(m1 == m2);
		}
	
		//
		// Binary Operators
		//
		friend Matrix<T> operator* ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
#if NSL_RANGE_CHECK
			if ( m1.cols() != m2.rows()) 
			{
				throw gsl_error( "Matrix<T> operator*: Matricies not compatible for multiply.",
					__FILE__, __LINE__, GSL_EINVAL);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize(m1.rows(), m2.cols());
			T result;
			for ( size_t i = 0; i < m1.rows(); ++i)
			{
				for ( size_t j = 0; j < m2.cols(); ++j) 
				{
					result = 0;
					for ( size_t k = 0; k < m1.cols(); ++k) 
					{
						
						result += m1.cref(i,k) * m2.cref(k,j);
					}
					cacheMatrix0_(i,j) = result;
				}
			}
			return cacheMatrix0_;
		};
		
		friend Matrix<T> operator* ( const GneMatrix<T> &m1, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ *= t;
			return cacheMatrix0_;
		};
		friend Matrix<T> operator* ( const T &t, const GneMatrix<T> &m1)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ *= t;
			return cacheMatrix0_;
		}
		
		friend valarray<T> operator* ( const GneMatrix<T> &m, const valarray<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( m.cols() != v.size())
			{
				throw gsl_error( "Matrix<T>::operator*: Matrix and Vector not compatible.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			valarray<T> vTemp( m.rows());
			T result;
			for ( size_t i = 0; i < m.rows(); ++i) {
				result = T(0);
				for ( size_t j = 0; j < m.cols(); ++j) {
					result += m(i,j) * v[j];
				}
				vTemp[i] = result;
			}
			return vTemp;
		}
		
		friend valarray<T> operator* ( const valarray<T> &v, const GneMatrix<T> &m)
		{
#if NSL_RANGE_CHECK
			if ( m.rows() != v.size()) 
			{
				throw gsl_error( "Matrix<T>::operator*: Vector and Matrix not compatible.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			valarray<T> vTemp( m.cols());
			T result;
			for ( size_t j = 0; j < m.cols(); ++j)
			{
				result = T(0);
				for ( size_t i = 0; i < m.rows(); ++i) 
				{
					result += v[i] * m(i,j);
				}
				vTemp[j] = result;
			}
			return vTemp;
		}

		friend Matrix<T> operator/ ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
			return m1 * !m2;
		};
		
		friend Matrix<T> operator/ ( const GneMatrix<T> &m1, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ /= t;
			return cacheMatrix0_;
		};
		
		friend Matrix<T> operator/ ( const T &t, const GneMatrix<T> &m1)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_.inverse();
			cacheMatrix0_ *= t;
			return cacheMatrix0_;
		}
		
		friend valarray<T> operator/ ( const GneMatrix<T> &m, const valarray<T> &v)
		{
			return m * (T(1)/v);
		}
		
		friend valarray<T> operator/ ( const valarray<T> &v, const GneMatrix<T> &m)
		{
			return v * !m;
		}

		friend Matrix<T> operator+ ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ += m2;
			return cacheMatrix0_;
		};
		
		friend Matrix<T> operator+ ( const GneMatrix<T> &m1, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ += t;
			return cacheMatrix0_;
		};
		friend Matrix<T> operator+ ( const T &t, const GneMatrix<T> &m1)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			ValMatrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ += t;
			return cacheMatrix0_;
		}

		friend Matrix<T> operator- ( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ -= m2;
			return cacheMatrix0_;
		};
		
		friend Matrix<T> operator- ( const GneMatrix<T> &m1, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ -= t;
			return cacheMatrix0_;
		};
		friend Matrix<T> operator- ( const T &t, const GneMatrix<T> &m1)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_  = m1;
			cacheMatrix0_ -= t;
			return cacheMatrix0_;
		}

		//
		// Misc. functions.
		//
// 		friend Matrix<T> abs( const GneMatrix<T> &m)
// 		{
// #ifdef NSL_SUPPORT_MULTITHREAD
// 			Matrix<T> cacheMatrix0_;
// #endif
// 			cacheMatrix0_.resize( m.rows(), m.cols());
// 			// for loop faster than using valarray abs() function.
// 			for (size_t i=0; i<rows(); i++)
// 			{
// 				for (size_t j=0; j<cols(); j++)
// 				{
// 					cacheMatrix0_(i,j) = nMath::abs(m.cref(i,j));
// 				}
// 			}
// 			return cacheMatrix0_;
// 		}
		
		friend Matrix<T> ceil( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());

			for (size_t i=0; i<rows(); i++)
			{
				for (size_t j=0; j<cols(); j++)
				{
					cacheMatrix0_(i,j) = ::ceil(m.cref(i,j));
				}
			}
			return cacheMatrix0_;
		}

		friend Matrix<T> floor( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());

			for (size_t i=0; i<rows(); i++)
			{
				for (size_t j=0; j<cols(); j++)
				{
					cacheMatrix0_(i,j) = ::floor(m.cref(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> pow( const GneMatrix<T> &m, size_t n)
		{
#if NSL_RANGE_CHECK
			if ( m.isSquare() == false) 
			{
				throw gsl_error( "pow(Matrix<T>): Matrix not square.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix1_;
			Matrix<T> cacheMatrix2_;
#endif
			
			if ( n == 0) 
			{
				cacheMatrix1_.resize( m.rows(), m.cols());
				cacheMatrix1_.unit();
				return cacheMatrix1_;
			}
			
			cacheMatrix2_  = m;
			bool firstTime = true;
			while( true) 
			{
				if ( n & 1) 
				{
					if ( firstTime) 
					{
						cacheMatrix1_ = cacheMatrix2_;
						firstTime     = false;
					}
					else 
					{
						cacheMatrix1_ *= cacheMatrix2_;
					}
				}
				n >>= 1;
				if ( n == 0) 
				{
					break;
				}
				cacheMatrix2_ *= cacheMatrix2_;
			}
			return cacheMatrix1_;
		}
		
		friend Matrix<T> acos( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray acos() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::acos(m(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> asin( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray asin() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::asin(m(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> atan( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray atan() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::atan(m(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> atan2( const GneMatrix<T> &m1, const GneMatrix<T> &m2)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m1.rows(), m1.cols());
			// for loop faster than using valarray atan2() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::atan2(m1(i,j), m2(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> atan2( const GneMatrix<T> &m, const T &v)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray atan2() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::atan2(m1(i,j), v);
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> atan2( const T &v, const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray atan2() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::atan2(m1(i,j), v);
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> cos( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray cos() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::cos(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> cosh( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray cosh() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::cosh(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}

		friend Matrix<T> exp( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray exp() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::exp(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> log( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray log() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::log(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> log10( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray log10() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::log10(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> sin( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray sin() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::sin(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> sinh( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray sinh() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::sinh(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> sqrt( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray sqrt() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::sqrt(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> tan( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray tan() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::tan(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}
		
		friend Matrix<T> tanh( const GneMatrix<T> &m)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_.resize( m.rows(), m.cols());
			// for loop faster than using valarray tanh() function.
			for (size_t i=0; i<m.rows(); i++)
			{
				for (size_t j=0; j<m.cols(); j++)
				{
					cacheMatrix0_(i,j) = ::tanh(m1(i,j));
				}
			}
			return cacheMatrix0_;
		}

		
		friend Matrix<T> concatenate( const GneMatrix<T> &m1, const GneMatrix<T> &m2 )
		{
			return concatenate( m1, m2, CC_HORIZONTAL );
		}
		
		//矩阵拼接，vertical or horizontal
		friend Matrix<T> concatenate( const GneMatrix<T> &m1, const GneMatrix<T> &m2,
			const Cc_Dir cc_dir )
		{
#if NSL_RANGE_CHECK
			if ( cc_dir == CC_HORIZONTAL && m1.rows() != m2.rows() ) 
			{
				throw gsl_error( "concatenate(m1,m2,CC_HORIZONTAL): Matrices must have same row size.",
					__FILE__, __LINE__, GSL_EBADLEN);
			} 
			else if ( cc_dir == CC_VERTICAL && m1.cols() != m2.cols() ) 
			{
				throw gsl_error( "concatenate(m1,m2,CC_VERTICAL): Matrices must have same column size.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			if ( cc_dir == CC_HORIZONTAL ) 
			{
				cacheMatrix0_.resize( m1.rows(), m1.cols() + m2.cols());
				for ( size_t i = 0; i < m1.rows(); ++i)
				{
					{
						for ( size_t j = 0; j < m1.cols(); ++j) 
						{
							cacheMatrix0_(i,j) = m1.cref(i,j);
						}
					}
					{
						for ( size_t j = 0; j < m2.cols(); ++j) 
						{
							cacheMatrix0_(i,j+m1.cols()) = m2.cref(i,j);
						}
					}
				}
			} else {
				cacheMatrix0_.resize( m1.rows() + m2.rows(), m1.cols());
				for ( size_t j = 0; j < m1.cols(); ++j)
				{
					{
						for ( size_t i = 0; i < m1.rows(); ++i) 
						{
							cacheMatrix0_(i,j) = m1.cref(i,j);
						}
					}
					{
						for ( size_t i = 0; i < m2.rows(); ++i) 
						{
							cacheMatrix0_(i+m1.rows(),j) = m2.cref(i,j);
						}
					}
				}
			}
			return cacheMatrix0_;
		}

		friend void swap( AlmostMatrix<T> &m1, AlmostMatrix<T> &m2)
		{
#if NSL_RANGE_CHECK
			if (m1.rows() != m2.rows() || m1.cols()!=m2.cols())
			{
				throw gsl_error("swap: Matrix size not comformat", __FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_(m1.rows(), m1.cols());
			assign(m1, cacheMatrix0_);
			assign(m2, m1);
			assign(cacheMatrix0_, m2);
#else
			cacheMatrix0_.resize(m1.rows(), m1.cols());
			assign(m1, cacheMatrix0_);
			assign(m2, m1);
			assign(cacheMatrix0_, m2);
#endif
		}
	};
	
	template <typename T>
	class ConstMatrixView : public GneMatrix<T>  //常视图
	{
	private:
		const valarray<T>* pMatrixData_;
		mslice viewSlice_;
		mslice orgSlice_;
		
	public:
		
		
		//重写虚函数
		virtual const mslice& getMatrixSlice()	const {return viewSlice_ ;}
		virtual const valarray<T>& constVal()	const {return *pMatrixData_ ;}
		virtual const mslice& getOrgSlice()		const {return orgSlice_ ;}
		
		//constructor
		ConstMatrixView()											:pMatrixData_(NULL), viewSlice_(0,0,0,0), orgSlice_(viewSlice_)
		{}
		
		virtual ~ConstMatrixView()									
		{ 
			const_cast<const valarray<T>*>(pMatrixData_) = NULL;
		}
		
		ConstMatrixView(const valarray<T> & val)					:pMatrixData_(&val), viewSlice_(0, 0, 1, val.size()), orgSlice_(viewSlice_)
		{}

  		ConstMatrixView(const valarray<T> & val, size_t i, size_t j)	:pMatrixData_(&val), viewSlice_(0,0, i, j), orgSlice_(viewSlice_)
  		{}
		

		//下面两个构造函数是为了
		//Matrix			 ---> ConstMatrixView
		//MatrixView		 ---> ConstMatrixView
		//ConstMatrixView	 ---> ConstMatrixView
		ConstMatrixView(const GneMatrix<T> & gmtx)		:pMatrixData_(&gmtx.constVal()), viewSlice_(gmtx.getMatrixSlice()), orgSlice_(gmtx.getOrgSlice())
		{}
		
		ConstMatrixView(const GneMatrix<T> & gmtx, const mslice& s)	:pMatrixData_(&gmtx.constVal()), viewSlice_(s), orgSlice_(gmtx.getOrgSlice())
		{}

		//subMatrix , 子矩阵常视图
		ConstMatrixView<T> subMatrix	(const mslice & s)		const 
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("subMatrix(mslice): mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}

		ConstMatrixView<T> operator []	(const mslice & s)		const
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}
		ConstMatrixView<T> operator ()	(const mslice & s)		const 
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}

		// 返回第i行的常视图
		ConstVectorView<T> operator [](size_t i) const { return ConstVectorView<T>( constVal(), slice( index(i, 0) , viewSlice_.cols(), 1 ));}
		ConstVectorView<T> operator ()(size_t i) const { return ConstVectorView<T>( constVal(), slice( index(i, 0) , viewSlice_.cols(), 1 ));}
		
		T operator ()(size_t i, size_t j) const { return constVal()[index(i, j)] ;}
		
		ConstVectorView<T> diagView( int i = 0)	 const 
		{
#if	NSL_RANGE_CHECK
			if ( i >= cols() || rows() + i < 1) 
			{
				throw gsl_error( "Matrix<T>::diag(i): Index i out of range.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			size_t start;
			size_t length;
			if ( i >=0) 
			{
				start  = i;
				length = cols() - i;
				if ( length > rows()) 
				{
					length = rows();
				}
			}
			else {
				start  = -i * cols();
				length = rows() + i;
				if ( length > cols()) 
				{
					length = cols();
				}
			}
			return ConstVectorView<T>(constVal(),slice( start, length, cols() + 1));
		}
	};
	
	template <typename T>
	class AlmostMatrix : public GneMatrix<T> //可变矩阵基类，继承于 GneMatrix, 被 MatrixView 和 Matrix 继承
	{
	
	public:
		//constructor and destructor
		AlmostMatrix()
		{}
		
		virtual ~AlmostMatrix()
		{}
		
		using GneMatrix<T>::getMatrixSlice;
		using GneMatrix<T>::getOrgSlice;
		
		//override this function in MatrixView and Matrix
		virtual valarray<T>& val() = 0;
		virtual gDecompositionInfo<T>& getGInfo() = 0;
		//MatrixView and Matrix 的共有函数

		//subMatrix , 子矩阵常视图
		ConstMatrixView<T> subMatrix	(const mslice & s)		const
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}
		ConstMatrixView<T> operator []	(const mslice & s)		const
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}

		ConstMatrixView<T> operator ()	(const mslice & s)		const
		{
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return ConstMatrixView<T>(*this, s) ;
		}

		//subMatrix , 子矩阵可变视图
		MatrixView<T> subMatrix	 ( const mslice & s)		
		{ 
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return MatrixView<T>(*this, s) ;
		}	

		MatrixView<T> operator[] ( const mslice & s)		
		{ 
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return MatrixView<T>(*this, s) ;
		}

		MatrixView<T> operator() ( const mslice & s)		
		{ 
#if NSL_RANGE_CHECK
			if (!hasSubMatrix(s))
			{
				throw gsl_error("operator[mslice]: mslice out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return MatrixView<T>(*this, s) ;
		}
		
		// 返回第i行的常视图
		ConstVectorView<T> operator [](size_t i) const { return ConstVectorView<T>( constVal(), slice( index(i, 0) , getMatrixSlice().cols(), 1 ));}
		ConstVectorView<T> operator ()(size_t i) const { return ConstVectorView<T>( constVal(), slice( index(i, 0) , getMatrixSlice().cols(), 1 ));}
		
		// 返回第i行的可变视图
		VectorView<T> operator [](size_t i)  { return VectorView<T>( val(), slice( index(i, 0) , getMatrixSlice().cols(), 1 ));}
		VectorView<T> operator ()(size_t i)  { return VectorView<T>( val(), slice( index(i, 0) , getMatrixSlice().cols(), 1 ));}
		
		//access function
		T& operator ()(size_t i, size_t j)  	  { return val()[index(i, j)] ;}
		T  operator ()(size_t i, size_t j)	const { return constVal()[index(i, j)] ;}

		inline T& ref(size_t i, size_t j)  		  { return val()[index(i, j)] ;}
		
		

		//
		// The default methods.
		//
		valarray<T> solve( const valarray<T> &v) const { return gSolve( v);}
		inline T det()							 const { return gDeterminant();}
		T        rank()							 const { return gRank();}
		T        condition()						 const { return svCondition();}

		//随机矩阵
		void rand( T dRandMin = T(-1), T dRandMax = T(1), int randSeed = 0)
		{
			T dR = dRandMax - dRandMin;
			if ( randSeed == 0) 
			{
				randSeed = (unsigned) time(0);
			}
			srand( randSeed);
			for ( size_t i = 0; i < rows(); ++i) 
			{
				for (size_t j=0; j<cols(); j++)
				{
					ref(i,j) = ( ::rand() * dR / RAND_MAX) + dRandMin;
				}
			}
		}

		VectorView<T> diagView( int i = 0)	
		{
#if	NSL_RANGE_CHECK
			if ( i >= cols() || rows() + i < 1) 
			{
				throw gsl_error( "Matrix<T>::diag(i): Index i out of range.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			size_t start;
			size_t length;
			if ( i >=0) 
			{
				start  = i;
				length = cols() - i;
				if ( length > rows()) 
				{
					length = rows();
				}
			}
			else 
			{
				start  = -i * cols();
				length = rows() + i;
				if ( length > cols()) 
				{
					length = cols();
				}
			}
			
			return VectorView<T>(val(),slice( start, length, cols() + 1));
		}

		ConstVectorView<T> diagView( int i = 0) const
		{
#if	NSL_RANGE_CHECK
			if ( i >= cols() || rows() + i < 1) 
			{
				throw gsl_error( "Matrix<T>::diag(i): Index i out of range.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			size_t start;
			size_t length;
			if ( i >=0) 
			{
				start  = i;
				length = cols() - i;
				if ( length > rows()) 
				{
					length = rows();
				}
			}
			else {
				start  = -i * cols();
				length = rows() + i;
				if ( length > cols()) 
				{
					length = cols();
				}
			}
			return ConstVectorView<T>(constVal(),slice( start, length, cols() + 1));
		}
		void unit()   
		{ 
			zero(); 
			diagView(0) = T(1) ;
		}
		void zero()  
		{ 
			for (size_t i=0; i<rows(); i++)
			{
				for (size_t j=0; j<cols(); j++)
				{
					ref(i,j) = T(0);
				}
			}
		
		}
		void null()               			
		{ 
			zero() ;
		}
		
		void copy(const GneMatrix<T> &amtx)
		{
			for(size_t i=0; i<rows(); i++)
			{
				for (size_t j=0; j<cols(); j++)
				{
					ref(i, j) = amtx.cref(i, j);
				}
				
			}
		}
		
		void copy(const valarray<T> & val)
		{
			size_t _size = size();
			size_t i=0;
			for(i=0; i<_size; i++)
			{
				ref(i, j) = val(index(i, j));
			}
		}
		
		void copy(const vector<T> & vec)
		{
			size_t _size = size();
			size_t i=0;
			for(i=0; i<_size; i++)
			{
				ref(i, j) = vec(index(i, j));
			}
		}
		
		//矩阵求逆
		void inverse()
		{
#if NSL_RANGE_CHECK
			if (!isSquare())
				throw gsl_error( "Matrix<T>::inverse(): Not a square matrix.",
				__FILE__, __LINE__, GSL_EINVAL);
#endif
			
			size_t i, j, k;
			double currentMax, currentValue;
			T scale(0);
			
			valarray<size_t> rowIndex( rows());
			for ( i = 0; i < rows(); ++i)
				rowIndex[i] = i;
			for ( k = 0; k < rows(); ++k) {
				i = k;
				currentMax = nMath::abs( (*this)(k,k));
				for ( j = k + 1; j < rows(); ++j) {
					if ((currentValue = nMath::abs(ref(j,k))) > currentMax) 
					{
						currentMax = currentValue;
						i  = j;
					}
				}
#if NSL_RANGE_CHECK
				if ( currentMax < nMath::ZERO) 
				{
					throw gsl_error( "Matrix<T>::inverse(): Singular matrix.",
						__FILE__, __LINE__, GSL_EINVAL);
				}
#endif
				if ( i != k)
				{
					std::swap( rowIndex[k], rowIndex[i]);
					for ( j = 0; j < cols(); ++j) {
						std::swap( ref(k,j), ref(i,j));
					}
				}
				scale = T(1) / (*this)(k,k);
				(*this)(k,k) = T(1);
				for ( j = 0; j < cols(); ++j) {
					(*this)(k,j) *= scale;
				}
				for ( i = 0; i < rows(); ++i) {
					if ( i != k) {
						scale        = (*this)(i,k);
						(*this)(i,k) = T(0);
						for ( j = 0; j < cols(); ++j) {
							(*this)(i,j) -= scale * (*this)(k,j);
						}
					}
				}
			}
			
			for ( j = 0; j < rows(); ++j) {
				if ( j != rowIndex[j]) {
					k = j + 1;
					while ( j != rowIndex[k]) {
						++k;
					}
					for ( i = 0; i < rows(); ++i) {
						std::swap( (*this)(i,j), (*this)(i,k));
					}
					std::swap( rowIndex[j], rowIndex[k]);
				}
			}
		}

		gDecompositionInfo<T> gDecomposition( valarray<size_t> &rowIndex)
		{
#if NSL_RANGE_CHECK
			if ( rowIndex.size() != rows())
			{
				throw gsl_error( "Matrix<T>::gDecomposition(): Index vector not compatible with matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			size_t i, ii, j;
			size_t rowI, rowII, iMax;
			T scale, dValue, dMax;
			getGInfo().rank       = 0;
			getGInfo().determinant = T(1);
			for ( i = 0; i < rows(); ++i) 
			{
				rowIndex[i] = i;
			}
			for ( i = 0; i < rows() - 1; ++i) 
			{
				// Find maximum value in column from this row down.
				iMax = i;
				rowI = rowIndex[iMax];
				dMax = nMath::abs( (*this)( rowI, i));
				for ( ii = i + 1; ii < rows(); ++ii) 
				{
					rowII  = rowIndex[ii];
					dValue = nMath::abs( (*this)( rowII, i));
					if (  dValue > dMax) {
						dMax = dValue;
						iMax = ii;
					}
				}
				// Swap rows by index;
				if ( iMax != i)
				{
					ii                = rowIndex[iMax];
					rowIndex[iMax]    = rowIndex[i];
					rowIndex[i]       = ii;
					getGInfo().determinant = -getGInfo().determinant;
				}
				// Scale max value row.
				rowI   = rowIndex[i];
				dValue = (*this)( rowI, i);
				if ( nMath::abs( dValue) < nMath::ZERO)
				{
					// Singular Matrix.
					return gDecompositionInfo<T>( i, T(0));
				}
				getGInfo().determinant *= dValue;
				scale = T(1) / dValue;
				for ( j = i; j < cols(); ++j)
				{
					(*this)( rowI, j) *= scale;
				}
				// Pivot other rows.
				for ( ii = i + 1; ii < rows(); ++ii) 
				{
					rowII = rowIndex[ii];
					scale = (*this)( rowII, i);
					for ( j = i + 1; j < cols(); ++j) 
					{
						(*this)( rowII, j) -= (*this)( rowI, j) * scale;
					}
				}
			}
			dValue = (*this)( rowIndex[i], i);
			if ( nMath::abs( dValue) < nMath::ZERO)
			{
				// Singular Matrix.
				return gDecompositionInfo<T>( i, T(0));
			}
			getGInfo().rank         = nMath::min( rows(), cols());
			getGInfo().determinant *= dValue;
			return getGInfo();
		}

		valarray<T> gVectorBackSubstitution( const valarray<size_t> &rowIndex) const
		{
			// Back substitute and reorder for solution vector.
			size_t i, j, rowI;
			valarray<T> vr( rows());
			// Solve for last variable.
			i     = rows() - 1;
			rowI  = rowIndex[i];
			vr[i] = (*this)( rowI, rows()) / (*this)( rowI, i);
			// Solve for all other variables.
			while ( i--)
			{
				rowI  = rowIndex[i];
				vr[i] = (*this)( rowI, rows());
				for ( j = i + 1; j < rows(); ++j) {
					vr[i] -= (*this)( rowI, j) * vr[j];
				}
			}
			return vr;
		}

		Matrix<T> gMatrixBackSubstitution( const valarray<size_t> &rowIndex) const
		{
			// Back substitute and reorder for solution matrix.
			size_t i, j, jj, rowI;
			T dValue;
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix1_;
#endif
			cacheMatrix1_.resize( rows(), cols() - rows());
			// Solve for last variable.
			i      = rows() - 1;
			rowI   = rowIndex[i];
			dValue = T(1) / (*this)( rowI, i);
			for ( jj = 0; jj < cacheMatrix1_.cols(); ++jj) {
				cacheMatrix1_( i, jj) = (*this)( rowI, rows() + jj) * dValue;
			}
			// Solve for all other variables.
			while ( i--) 
			{
				rowI  = rowIndex[i];
				for ( jj = 0; jj < cacheMatrix1_.cols(); ++jj) 
				{
					cacheMatrix1_( i, jj) = (*this)( rowI, rows() + jj);
				}
				for ( j = i + 1; j < rows(); ++j) 
				{
					dValue = (*this)( rowI, j);
					for ( jj = 0; jj < cacheMatrix1_.cols(); ++jj) 
					{
						cacheMatrix1_( i, jj) -= dValue * cacheMatrix1_( j, jj);
					}
				}
			}
			return cacheMatrix1_;
		}

		void gBuildDecompositionMatrix( const valarray<T> &v, Matrix<T> &md) const
		{
			md.resize( rows(), cols() + 1);
			for ( size_t i = 0; i < rows(); ++i) {
				for ( size_t j = 0; j < cols(); ++j) {
					md( i, j) = (*this)( i, j);
				}
				md( i, cols()) = v[i];
			}
		}

		void gBuildDecompositionMatrix( const GneMatrix<T> &m, Matrix<T> &md) const
		{
			md.resize( rows(), cols() + m.cols());
			for ( size_t i = 0; i < rows(); ++i)
			{
				size_t j;
				for (j = 0; j < cols(); ++j) 
				{
					md( i, j) = (*this)( i, j);
				}
				for (j = 0; j < m.cols(); ++j) 
				{
					md( i, cols() + j) = m.cref( i, j);
				}
			}
		}

		valarray<T> gSolve( const valarray<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false) 
			{
				throw gsl_error( "Matrix<T>::gSolve(): Not a square matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
			if ( v.size() != rows())
			{
				throw gsl_error( "Matrix<T>::gSolve(): Vector not compatible with matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			valarray<size_t> rowIndex( rows());
			gBuildDecompositionMatrix( v, cacheMatrix0_);
			getGInfo() = cacheMatrix0_.gDecomposition( rowIndex);
			if ( getGInfo().rank == nMath::min( rows(), cols())) {
				return cacheMatrix0_.gVectorBackSubstitution( rowIndex);
			}
			return valarray<T>( T(0), rows());
		}

		Matrix<T> gSolve( const GneMatrix<T> &m)
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false) 
			{
				throw gsl_error( "Matrix<T>::gSolve(): Not a square matrix.",
						__FILE__, __LINE__, GSL_EBADLEN);
			}
			if ( m.rows() != rows()) 
			{
				throw gsl_error( "Matrix<T>::gSolve(): Input matrix not compatible.",
						__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			valarray<size_t> rowIndex( rows());
			gBuildDecompositionMatrix( m, cacheMatrix0_);
			getGInfo() = cacheMatrix0_.gDecomposition( rowIndex);
			if ( getGInfo().rank == nMath::min( rows(), cols())) 
			{
				return cacheMatrix0_.gMatrixBackSubstitution( rowIndex);
			}
			cacheMatrix0_.resize( m.rows(), m.cols());
			cacheMatrix0_ = T(0);
			return cacheMatrix0_;
		}

		T gDeterminant() const
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false) 
			{
				throw gsl_error( "Matrix<T>::gDeterminant(): Not a square matrix.",
						__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = (*this);
			valarray<size_t> rowIndex( rows());
			return cacheMatrix0_.gDecomposition( rowIndex).determinant;
		}
		
		T gRank() const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = (*this);
			valarray<size_t> rowIndex(  rows());
			return cacheMatrix0_.gDecomposition( rowIndex).rank;
		}

		//
		// Matrix LU functions.
		//
		T luDecomposition( valarray<size_t> &rowIndex)
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false) 
			{
				throw gsl_error( "Matrix<T>::luDecomposition(): Not a square matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			
			size_t i, j, k;
			double currentMax, currentValue;
			T scale;
			
			T determinant = T(1);
			if ( rows() != rowIndex.size()) 
			{
				rowIndex.resize( rows());
			}
			for ( i = 0; i < rows(); ++i)
			{
				rowIndex[i] = i;
			}
			for ( k = 0; k < rows() - 1; ++k) 
			{
				j = k;
				currentMax = nMath::abs( (*this)(rowIndex[k],k));
				for ( i = k + 1; i < rows(); ++i) 
				{
					if ((currentValue = nMath::abs( (*this)(rowIndex[i],k))) > currentMax) 
					{
						currentMax = currentValue;
						j = i;
					}
				}
				if ( j != k) {
					std::swap( rowIndex[j], rowIndex[k]);
					determinant = -determinant;
				}
				if ( nMath::abs( (*this)(rowIndex[k],k)) < nMath::ZERO) 
				{
					return T(0);
				}
				determinant *= (*this)(rowIndex[k],k);
				for ( i = k + 1; i < rows(); ++i) 
				{
					scale = (*this)(rowIndex[i],k) /= (*this)(rowIndex[k],k);
					for ( j = k + 1; j < rows(); ++j) 
					{
						(*this)(rowIndex[i],j) -= scale * (*this)(rowIndex[k],j);
					}
				}
			}
			determinant *= (*this)(rowIndex[k],k);
			return determinant;
		}

		valarray<T> luBackSubstitution( const valarray<size_t> &rowIndex, const valarray<T> &v) const
		{
			size_t i, j, k, ri;
			bool equalToZero = true;
			
#if NSL_RANGE_CHECK
			if ( v.size() != rows_) 
			{
				throw gsl_error( "VMatrix<T>::luBackSubstitution(): Incorrect valarray size.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
			if ( rowIndex.size() != rows_)
			{
				throw gsl_error( "Matrix<T>::luBackSubstitution(): Incorrect rowIndex size.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			
			valarray<T> solution( rows());
			k = 0;
			for ( i = 0; i < rows(); ++i) 
			{
				ri = rowIndex[i];
				solution[i] = v[ri];
				if ( !equalToZero) 
				{
					for ( j = k; j < i; ++j)
					{
						solution[i] -= (*this)(ri,j) * solution[j];
					}
				}
				else if ( vm_math::abs( solution[i]) > vm_math::ZERO) 
				{
					k = i;
					equalToZero = false;
				}
			}
			for ( i = rows() - 1; ; --i) 
			{
				ri = rowIndex[i];
				for ( j = i + 1; j < rows(); ++j) 
				{
					solution[i] -= (*this)(ri,j) * solution[j];
				}
				solution[i] /= (*this)(ri,i);
				if ( i == 0) break;
			}
			return solution;
		}

		valarray<T> luSolve( const valarray<T> &v) const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_= *this;
			valarray<size_t> rowIndex( rows());
#if NSL_RANGE_CHECK
			T determinant = cacheMatrix0_.luDecomposition( rowIndex);
			if ( nMath::abs( determinant) < nMath::ZERO) 
			{
				throw gsl_error( "Matrix<T>::luDecomposition(): Singular matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#else
			cacheMatrix0_.luDecomposition( rowIndex);
#endif
			return cacheMatrix0_.luBackSubstitution( rowIndex, v);
		}
		
		T luDeterminant() const
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false)
				throw gsl_error( "ValMatrix<T>::luDeterminant(): Not a square matrix.",
					__FILE__, __LINE__, GSL_EBADLEN);
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix0_;
#endif
			cacheMatrix0_ = *this;
			valarray<size_t> rowIndex( rows());
			return cacheMatrix0_.luDecomposition( rowIndex);
		}

			/** Singular Value Decomposition.
     * This function computes the singular value decomposition of a matrix 'a':
     *  [a] = [u][diag(w)][v]^T
     *
     * Base Matrix  a(m,n) : Matrix to decompose.
     * Base Matrix  u(m,n) : LH matrix of decomposition (left orthogonal
     *                       transformation) is returned to base matrix.
     * @return      v(n,n) : RH matrix of decomposition (transpose of right
     *                       orthogonal transformation (**not v^T **)).
     * @return      w[n]   : Diagonal singular values vector.
     * @code a.svDecomposition( w, v); @endcode
     */
		void svDecomposition( AlmostMatrix<T> &v, valarray<T> &w, T tolerance = nMath::ZERO)
		{
			if ( v.rows() != cols() || v.cols() != cols()) 
			{
				v.resize( cols(), cols());
			}
			if ( w.size() != cols()) 
			{
				w.resize( cols());
			}
			
			valarray<T> rv1( cols());
			
			size_t flag;
			size_t i, j, jj, k, l(0), nm;
			size_t iterations;
			T c, f, h, s, x, y, z, tmp;
			T g(0), scale(0), anorm(0);
			for ( i = 0; i < cols(); ++i) 
			{
				l      = i + 1;
				rv1[i] = scale * g;
				g      = T(0);
				s      = T(0);
				scale  = T(0);
				if ( i < rows()) 
				{
					for ( k = i; k < rows(); ++k) 
					{
						scale += nMath::abs( (*this)(k,i));
					}
					if ( scale > tolerance) 
					{
						for ( k = i; k < rows(); ++k) 
						{
							tmp = (*this)(k,i) /= scale;
							s  += tmp * tmp;
						}
						f = (*this)(i,i);
						g = -dstomath::sign( ::sqrt(s), f);
						h = f * g - s;
						(*this)(i,i) = f - g;
						
						for ( j = l; j < cols(); ++j) 
						{
							s = T(0);
							for ( k = i; k < rows(); ++k) 
							{
								s += (*this)(k,i) * (*this)(k,j);
							}
							f = s / h;
							for ( k = i; k < rows(); ++k) 
							{
								(*this)(k,j) += f * (*this)(k,i);
							}
						}
						for ( k= i; k < rows(); ++k) 
						{
							(*this)(k,i) *= scale;
						}
					}
				}
				w[i]   = scale * g;
				g      = T(0);
				s      = T(0);
				scale  = T(0);
				if ( i < rows() && i != cols()-1) 
				{
					for ( k = l; k < cols(); ++k)
					{
						scale += nMath::abs((*this)(i,k));
					}
					if ( scale > tolerance) 
					{
						for ( k = l; k < cols(); ++k) 
						{
							tmp = (*this)(i,k) /= scale;
							s  += tmp * tmp;
						}
						f = (*this)(i,l);
						g = -dstomath::sign( ::sqrt(s), f);
						h = f * g - s;
						(*this)(i,l) = f - g;
						for ( k = l; k < cols(); ++k) 
						{
							rv1[k] = (*this)(i,k) / h;
						}
						for ( j = l; j < rows(); ++j) 
						{
							s = T(0);
							for ( k = l; k < cols(); ++k)
							{
								s += (*this)(j,k) * (*this)(i,k);
							}
							for ( k = l; k < cols(); ++k) 
							{
								(*this)(j,k) += s * rv1[k];
							}
						}
						for ( k = l; k < cols(); ++k) 
						{
							(*this)(i,k) *= scale;
						}
					}
				}
				anorm = nMath::max( anorm, nMath::abs( w[i]) + nMath::abs( rv1[i]));
			}
			
			i = cols();
			do 
			{
				--i;
				if ( i < (cols() - 1)) 
				{
					if ( nMath::abs(g) > tolerance) 
					{
						for ( j = l; j < cols(); ++j) 
						{
							v(j,i) = ((*this)(i,j) / (*this)(i,l)) / g;
						}
						for ( j = l; j < cols(); ++j) 
						{
							s = T(0);
							for ( k = l; k < cols(); ++k)
							{
								s += (*this)(i,k) * v(k,j);
							}
							for ( k = l; k < cols(); ++k) 
							{
								v(k,j) += s * v(k,i);
							}
						}
					}
					for ( j = l; j < cols(); ++j) 
					{
						v(i,j) = v(j,i) = T(0);
					}
				}
				v(i,i) = T(1);
				g = rv1[i];
				l = i;
			} while ( i != 0);
			
			i = nMath::min( rows(), cols());
			do 
			{
				--i;
				l = i + 1;
				g = w[i];
				for ( j = l; j < cols(); ++j) 
				{
					(*this)(i,j) = T(0);
				}
				if ( nMath::abs(g) > tolerance) 
				{
					g = T(1) / g;
					for ( j = l; j < cols(); ++j)
					{
						s = T(0);
						for ( k = l; k < rows(); ++k)
						{
							s += (*this)(k,i) * (*this)(k,j);
						}
						f = ( s / (*this)(i,i)) * g;
						for ( k = i; k < rows(); ++k)
						{
							(*this)(k,j) += f * (*this)(k,i);
						}
					}
					for (j = i; j < rows(); ++j) 
					{
						(*this)(j,i) *= g;
					}
				}
				else 
				{
					for ( j = i; j < rows(); ++j) 
					{
						(*this)(j,i) = T(0);
					}
				}
				++(*this)(i,i);
			} while (i != 0);
			
			k = cols();
			do 
			{
				--k;
				for ( iterations = 1; iterations <= 30; ++iterations) 
				{
					flag = 1;
					l = k + 1;
					do 
					{
						--l;
						nm = l - 1;
						if ( nMath::abs( rv1[l]) < tolerance) 
						{
							flag = 0;
							break;
						}
						if ( nMath::abs( w[nm]) < tolerance)
						{
							break;
						}
					} while ( l != 0);
					if ( flag) 
					{
						c = T(0);
						s = T(1);
						for ( i = l; i <= k; ++i) 
						{
							f      = s * rv1[i];
							rv1[i] = c * rv1[i];
							if ( nMath::abs(f) < tolerance) 
							{
								break;
							}
							g    = w[i];
							h    = nMath::pythag( f, g);
							w[i] = h;
							h    = 1.0 / h;
							c    = g * h;
							s    = -f * h;
							for ( j = 0; j < rows(); ++j) 
							{
								y = (*this)(j,nm);
								z = (*this)(j,i);
								(*this)(j,nm) = y * c + z * s;
								(*this)(j,i)  = z * c - y * s;
							}
						}
					}
					z = w[k];
					if ( l == k) 
					{
						if ( z < T(0)) 
						{
							w[k] = -z;
							for ( j = 0; j < cols(); ++j) 
							{
								v(j,k) = -v(j,k);
							}
						}
						break;
					}
					
					if ( iterations == 30)
					{
#if NSL_RANGE_CHECK
						throw gsl_error( "Matrix<T>::svDecomposition: No convergence after 30 iterations.",
							__FILE__, __LINE__, GSL_EINVAL);
#endif
						break;
					}
					
					x  = w[l];
					nm = k - 1;
					y  = w[nm];
					g  = rv1[nm];
					h  = rv1[k];
					f  = ((y-z) * (y+z) + (g-h) * (g+h)) / (2.0 * h * y);
					g  = nMath::pythag( f, T(1));
					f  = ((x-z) * (x+z) + h * ((y / (f+nMath::sign( g, f))) - h)) / x;
					c  = T(1);
					s  = T(1);
					for ( j = l; j <= nm; ++j) 
					{
						i      = j + 1;
						g      = rv1[i];
						y      = w[i];
						h      = s * g;
						g      = c * g;
						z      = nMath::pythag( f, h);
						rv1[j] = z;
						c      = f / z;
						s      = h / z;
						f      = x * c + g * s;
						g      = g * c - x * s;
						h      = y * s;
						y     *= c;
						for ( jj = 0; jj < cols(); ++jj) 
						{
							x = v(jj,j);
							z = v(jj,i);
							v(jj,j) = x * c + z * s;
							v(jj,i) = z * c - x * s;
						}
						z    = nMath::pythag( f, h);
						w[j] = z;
						if ( nMath::abs(z) > tolerance) 
						{
							z = 1.0 / z;
							c = f * z;
							s = h * z;
						}
						f = c * g + s * y;
						x = c * y - s * g;
						for ( jj = 0; jj < rows(); ++jj) 
						{
							y = (*this)(jj,j);
							z = (*this)(jj,i);
							(*this)(jj,j) = y * c + z * s;
							(*this)(jj,i) = z * c - y * s;
						}
					}
					rv1[l] = 0.0;
					rv1[k] = f;
					w[k]   = x;
        }
      } while ( k != 0);
      return;
		}

		valarray<T> svBackSubstitution( const AlmostMatrix<T> &v, const valarray<T> &w,
			const valarray<T> &b) const
		{
#if NSL_RANGE_CHECK
			if ( v.rows() != cols_ || v.cols() != cols_) 
			{
				throw gsl_error( "Matrix<T>::svBackSubstitution(): Incorrect Matrix size for v.".
					__FILE__, __LINE__, GSL_EINVAL);
			}
			if ( w.size() != cols_) 
			{
				throw gsl_error( "Matrix<T>::svBackSubstitution(): Incorrect valarray size for w.".
					__FILE__, __LINE__, GSL_EINVAL);
			}
			if ( b.size() != rows_)
			{
				throw gsl_error( "Matrix<T>::svBackSubstitution(): Incorrect valarray size for b.".
					__FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			
			//
			// Only multiply matrix elements for w[j] > ZERO.
			// This increases matrix speed over a Vect2 = Matrix * Vect1
			// and prevents a divide by ZERO later.
			//
			size_t i, j;
			valarray<T> tmp( cols_);
			for ( j = 0; j < cols_; ++j) 
			{
				T result(0);
				if ( nMath::abs( w[j]) > nMath::ZERO) 
				{
					for ( i = 0; i < rows_; ++i) 
					{
						result += (*this)(i,j) * b[i];
					}
					result /= w[j];
				}
				tmp[j] = result;
			}
			
			return (v * tmp);
		}

		/** Solve Set of Equations using Singular Value Decomposition.
		* This functions solves A.X = B for a vector X, where A is represented by
		* a set of matrices determined from singular value decomposition of A.
		* The function evaluates:
		*    [s] = [v] * [diag(1/w)] * [u]^T * [b]
		*
		* Base Matrix a(m,n) : Matrix to decompose and left const.
		* @param      b[m]   : RHS to be solved for.
		* @return     s[m]   : Solution vector.
		* @code s = a.svBackSubstitution( v, w, b); @endcode
		*/
		valarray<T> svSolve( const valarray<T> &b) const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			ValMatrix<T> cacheMatrix1_;
			ValMatrix<T> cacheMatrix2_;
#endif
			cacheMatrix1_= *this;
			valarray<T> w;
			cacheMatrix1_.svDecomposition( cacheMatrix2_, w);
			/*
			* This is a very simple test for rank deficient matrices, setting
			* the factors of the deficient basis functions to zero.  Some
			* better algorithm may be required in future.
			*/
			T cutoff = w.max() * T(1.0e-5);
			valarray<bool> tCut = ( w - cutoff ) > T(0);
			w = w * T(0) + static_cast<valarray<double> >( w[ tCut ] );
			return cacheMatrix1_.svBackSubstitution( cacheMatrix2_, w, b);
		}
		
		/** Covariance of Singular Value Decomposition.
		* This function evaluates the covariance matrix for the solution of
		* A.X = B achieved using singular value decomposition.
		* Base Matrix v(n,n) : RH matrix of decomposition.
		* @param      w[n]   : Singular values.
		* @return covar(n,n) : The covariance matrix.
		* @code covar = v.svCovarience( w); @endcode
		*/
		Matrix<T> svCovariance( const valarray<T> &w)
		{
#if NSL_RANGE_CHECK
			if ( isSquare() == false) 
			{
				throw gsl_error( "Matrix<T>::svCovariance(): Matrix must be square.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
			if ( w.size() != cols()) 
			{
				throw gsl_error( "Matrix<T>::svCovariance(): Incorrect valarray size for w.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			
#ifdef NSL_SUPPORT_MULTITHREAD
			ValMatrix<T> cacheMatrix1_;
			ValMatrix<T> cacheMatrix2_;
#endif
			cacheMatrix1_.resize( rows(), cols());
			cacheMatrix2_.resize( rows(), cols());
			for ( size_t i = 0; i < rows_; i++) 
			{
				if ( nMath::abs( w[i]) > nMath::ZERO) 
				{
					cacheMatrix1_(i) = (*this)(i) / w[i];
				}
				else {
					cacheMatrix1_(i) = T(0);
				}
			}
			for ( size_t j = 0; j < rows(); ++j) 
			{
				for ( size_t k = j; k < rows(); ++k) 
				{
					cacheMatrix2_(j,k) = dot( cacheMatrix1_[j], cacheMatrix1_[k]);
					if ( k != j) {
						cacheMatrix2_(k,j) = cacheMatrix2_(j,k);
					}
				}
			}
			return cacheMatrix2_;
		}

		size_t svRank() const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix1_;
			Matrix<T> cacheMatrix2_;
#endif
			cacheMatrix1_= *this;
			valarray<T> w( cols_);
			cacheMatrix1_.svDecomposition( cacheMatrix2_, w);
			size_t count = 0;
			for ( size_t i = 0; i < cols_; ++i) 
			{
				if ( nMath::abs( w[i]) > nMath::ZERO)
				{
					++count;
				}
			}
			return count;
		}
		
		T svCondition() const
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Matrix<T> cacheMatrix1_;
			Matrix<T> cacheMatrix2_;
#endif
			cacheMatrix1_= *this;
			valarray<T> w;
			cacheMatrix1_.svDecomposition( cacheMatrix2_, w);
			return w.max()/w.min();
		}

		CommaAssignmentForMatrix<T> operator<< ( const T &t)
		{
			return CommaAssignmentForMatrix<T>( *this, t);
		}

		
	};
		
	
	template <typename T>
	class MatrixView : public AlmostMatrix<T> //可变视图类
	{
		private:
		valarray<T>* pMatrixData_;
		mslice viewSlice_;
		mslice orgSlice_;
		mutable gDecompositionInfo<T> gInfo_;	
	public:
		//constructions and destruction
		MatrixView()
		{}
		
		virtual ~MatrixView()
		{const_cast<valarray<T>*>(pMatrixData_) = NULL;}
		
		MatrixView(valarray<T> & val)					:pMatrixData_(&val), viewSlice_(0, 0, 1, val.size()), orgSlice_(viewSlice_)
		{}

		MatrixView(valarray<T> & val, size_t i, size_t j) :pMatrixData_(&val), viewSlice_(0, 0, i, j), orgSlice_(viewSlice_)
		{}
		

		//下面两个构造函数是为了
		//Matrix			 ---> MatrixView
		//MatrixView		 ---> MatrixView
		MatrixView(AlmostMatrix<T> & amtx)		:pMatrixData_(&amtx.val()), viewSlice_(amtx.getMatrixSlice()), orgSlice_(amtx.getOrgSlice())
		{}
		
		MatrixView(AlmostMatrix<T> & amtx, const mslice& s)	:pMatrixData_(&amtx.val()), viewSlice_(s), orgSlice_(amtx.getOrgSlice())
		{}
		
		
		//重写虚函数
		virtual const mslice& getMatrixSlice()	const {return viewSlice_ ;}
		virtual const valarray<T>& constVal()	const {return *pMatrixData_ ;}
		virtual const mslice& getOrgSlice()		const {return orgSlice_ ;}
		virtual valarray<T>& val() 	 				  {return *pMatrixData_ ;}
		virtual gDecompositionInfo<T>& getGInfo() {return gInfo_;};
		
		//
		// Assignment operators.
		//	
		
		MatrixView<T>& operator= ( const AlmostMatrix<T> &amtx)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size())
			{
				throw gsl_error( "void operator=: Matrix size not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			copy(amtx);
			return *this;
		}
		
		MatrixView<T>& operator= ( const valarray<T> & val)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			copy(val);
			return *this;
		}
		
		// This is for A[i] = vector.
		MatrixView<T>& operator= ( const std::vector<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator=: vector not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			copy(v);
			return *this;
		}

		
		MatrixView<T>& operator= ( const T t[])
		{
			size_t k= 0;
			for ( size_t i = 0; i < rows(); ++i) 
			{
				for (size_t j=0; j<cols(); ++j)
				{
					ref(i, j) = t[k++];
				}
			}
			return *this;
		}
		
		MatrixView<T>& operator= ( const T &t)
		{
			for ( size_t i = 0; i < rows(); ++i) 
			{
				for (size_t j=0; j<cols(); ++j)
				{
					ref(i, j) = t
				}
			}
			return *this;
		}

			//
			// Computed assignment operators.
			//
			MatrixView<T>& operator*= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( cols() != m.rows())
				{
					throw gsl_error("Matrix<T>::operator*=: Matricies not compatible for multiply.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
				Matrix<T> cacheMatrix0_;
#endif
				cacheMatrix0_.resize(rows(), cols());
				T result;
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < m.cols(); ++j) 
					{
						result = T(0);
						for ( size_t k = 0; k < m.rows(); ++k) 
						{
							result += cref(i,k) * m(k,j);
						}
						cacheMatrix0_(i,j) = result;
					}
				}
				*this = cacheMatrix0_;
				return *this;
			}
			
			MatrixView<T>& operator*= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) *= t;
					}
				}
				return *this;
			}
			
			MatrixView<T>& operator/= ( const Matrix<T> &m)
			{
				*this = *this * !m;
				return *this;
			}
			
			MatrixView<T>& operator/= ( const T &t)
			{
				*this *= 1.0 / t;
				return *this;
			}
			
			MatrixView<T>& operator+= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( (cols() != m.cols()) || (rows() != m.rows())) 
				{
					throw gsl_error( "Matrix<T>::operator+=: Matricies not compatible for addition.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) += m(i,j);
					}
				}
				return *this;
			}
			
			MatrixView<T>& operator+= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) += t;
					}
				}
				return *this;
			}
			
			MatrixView<T>& operator-= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( (cols() != m.cols()) || (rows() != m.rows())) 
				{
					throw gsl_error( "Matrix<T>::operator-=: Matricies not compatible for addition.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) -= m(i,j);
					}
				}
				return *this;
			}
			
			MatrixView<T>& operator-= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) -= t;
					}
				}
				return *this;
			}
		
	};
	
	
	template <typename T>
		class Matrix : public AlmostMatrix<T> //矩阵类
	{

		private:
			valarray<T> matrixData_;
			mslice		mslice_;
			mslice		orgSlice_;
			mutable gDecompositionInfo<T> gInfo_;	
		public:
			
			//重写的虚函数
			virtual const mslice& getMatrixSlice() 	const {return mslice_ ;}
			virtual const valarray<T>& constVal()  	const {return matrixData_ ;}
			virtual const mslice& getOrgSlice()		const {return orgSlice_ ;}
			virtual valarray<T>& val() 	 				  {return matrixData_ ;}
			virtual gDecompositionInfo<T>& getGInfo() {return gInfo_;};
			//constructions and destruction
			Matrix()										:matrixData_(0), mslice_(0,0,0,0),orgSlice_(mslice_)
			{}
			
			virtual ~Matrix()	{}
			
			Matrix(size_t i, size_t j)						:matrixData_(T(0), i*j), mslice_(0, 0, i, j), orgSlice_(mslice_)
			{}
			
			Matrix(size_t i, size_t j, const T & val)		:matrixData_(val, i*j), mslice_(0, 0, i, j), orgSlice_(mslice_)
			{}
			
			Matrix(const valarray<T> & val)					:matrixData_(val), mslice_(0, 0, 1, val.size()), orgSlice_(mslice_)
			{}

			Matrix(const valarray<T> & val, size_t i, size_t j)	:matrixData_(val[slice(0, i*j, 1)]), mslice_(0, 0,i, j), orgSlice_(mslice_)
			{}
			
			//由数组构造
			
			Matrix(const T * p, size_t i, size_t j) 			: matrixData_(p, i*j), mslice_(0, 0, i, j),orgSlice_(mslice_)
			{}

			//下面两个构造函数是为了
			//Matrix			 ---> Matrix
			//MatrixView		 ---> Matrix
			Matrix(const GneMatrix<T> & gmtx)				:matrixData_(gmtx.constVal()), mslice_(gmtx.getMatrixSlice()), orgSlice_(gmtx.getOrgSlice())
			{}
			
			Matrix(const GneMatrix<T> & gmtx, const mslice& ms)	:matrixData_(gmtx.constVal()), mslice_(ms), orgSlice_(mslice_)
			{}
			
			
			void resize( size_t i, size_t j )		{ matrixData_.resize(i*j); mslice_ = mslice(0, 0, i, j); orgSlice_= mslice_; }
	
			//
			// Assignment operators.
			//
			
			Matrix<T>& operator= ( const GneMatrix<T> &m)
			{
				if ( &m != this) 
				{
					if (size() != m.size())
					{
						resize(m.rows(), m.cols());
					}
					copy( m);
				}
				mslice_ = mslice(0,0,m.rows(), m.cols());
				orgSlice_ = mslice_;
				return *this;
			}
			
			Matrix<T>& operator= ( const valarray<T> &v)
			{
#if NSL_RANGE_CHECK
				if ( size() != v.size())
				{
					throw gsl_error( "Matrix<T>::operator=: valarray has incorrect size.",__FILE__,
						__LINE__, GSL_EINVAL);
				}
#endif
				matrixData_ = v;
				mslice_ = mslice(0,0, 1, v.size());
				orgSlice_ = mslice_;
				return *this;
			}
			
			Matrix<T>& operator= ( const T t[])
			{
				for ( size_t i = 0; i < size(); ++i) 
				{
					matrixData_[i] = t[i];
				}
				return *this;
			}
			
			Matrix<T>& operator= ( const T &t)
			{
				matrixData_ = t;
				return *this;
			}

			
			//
			// Computed assignment operators.
			//
			Matrix<T>& operator*= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( cols() != m.rows())
				{
					throw gsl_error("Matrix<T>::operator*=: Matricies not compatible for multiply.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
				Matrix<T> cacheMatrix0_;
#endif
				cacheMatrix0_.resize(rows(), cols());
				T result;
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < m.cols(); ++j) 
					{
						result = T(0);
						for ( size_t k = 0; k < m.rows(); ++k) 
						{
							result += cref(i,k) * m.cref(k,j);
						}
						cacheMatrix0_(i,j) = result;
					}
				}
				*this = cacheMatrix0_;
				return *this;
			}
			
			Matrix<T>& operator*= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) *= t;
					}
				}
				return *this;
			}
			
			Matrix<T>& operator/= ( const Matrix<T> &m)
			{
				*this = *this * !m;
				return *this;
			}
			
			Matrix<T>& operator/= ( const T &t)
			{
				*this *= 1.0 / t;
				return *this;
			}
			
			Matrix<T>& operator+= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( (cols() != m.cols()) || (rows() != m.rows())) 
				{
					throw gsl_error( "Matrix<T>::operator+=: Matricies not compatible for addition.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) += m.cref(i,j);
					}
				}
				return *this;
			}
			
			Matrix<T>& operator+= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) += t;
					}
				}
				return *this;
			}
			
			Matrix<T>& operator-= ( const GneMatrix<T> &m)
			{
#if NSL_RANGE_CHECK
				if ( (cols() != m.cols()) || (rows() != m.rows())) 
				{
					throw gsl_error( "Matrix<T>::operator-=: Matricies not compatible for addition.",
						__FILE__, __LINE__, GSL_EBADLEN);
				}
#endif
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) -= m.cref(i,j);
					}
				}
				return *this;
			}
			
			Matrix<T>& operator-= ( const T &t)
			{
				for ( size_t i = 0; i < rows(); ++i) 
				{
					for ( size_t j = 0; j < cols(); ++j) 
					{
						ref(i, j) -= t;
					}
				}
				return *this;
			}
	};

#ifndef NSL_SUPPORT_MULTITHREAD
	// 矩阵运算缓存.
	template <typename T> Matrix<T> GneMatrix<T>::cacheMatrix0_;
	template <typename T> Matrix<T> GneMatrix<T>::cacheMatrix1_;
	template <typename T> Matrix<T> GneMatrix<T>::cacheMatrix2_;
#endif

	template <typename T> inline void assign(const AlmostMatrix<T> &m1, AlmostMatrix<T> &m2)
	{
		for (size_t i=0; i<m1.rows(); i++)
		{
			for (size_t j=0; j<m1.cols(); j++)
			{
				m2(i,j) = m1(i,j);
			}
				
		}
	}
	

	} // end of namespace gslcpp
#endif // NSL_MATRIX_H__