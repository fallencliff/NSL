
/********************************************************************
filename: 	VTest.h
author:		hu zhijian
created:	23:5:2010   19:39
brief:	
*********************************************************************/

#ifndef NSL_VECTOR_H__
#define NSL_VECTOR_H__


#include <NSL.h>
#include <valarray>
#include <cmath>
#include <float.h>
#include <vector>
#include <gsl_errno.h>
using std::valarray;
using std::slice;
using std::slice_array;
using std::vector;

namespace gslcpp
{
	
	
	namespace nMath
	{
		
		template <typename T> inline T abs( const T &t)               {return ((t < T(0)) ? -t : t);}
		template <typename T> inline T min( const T &t1, const T &t2) {return ((t1 < t2) ? t1 : t2);}
		template <typename T> inline T max( const T &t1, const T &t2) {return ((t1 > t2) ? t1 : t2);}
		const double ZERO = DBL_EPSILON * 100.0;
		
		//
		// General math functions.
		//
		template <typename T> inline T sign( const T &a, const T &b)
		{
			return ( b >= T(0) ? vm_math::abs(a) : -vm_math::abs(a));
		}
		
		template <typename T> inline T pythag( T a, T b)
		{
			T c;
			a = nMath::abs( a);
			b = nMath::abs( b);
			if ( a > b) 
			{
				c = b / a;
				return ( a * ::sqrt( T(1) + c * c));
			}
			if ( b < nMath::ZERO) 
			{
				return T(0);
			}
			c = a / b;
			return ( b * ::sqrt( T(1) + c * c));
		}
	}
	
	//
	// vector to valarray conversion.
	//
	template <typename T> inline valarray<T> getValArrayFromVector( const std::vector<T> &v)
	{
		valarray<T> vr( v.size());
		for ( size_t i = 0; i < v.size(); ++i)
		{
			vr[i] = v[i];
		}
		return vr;
	}
	
	//前向声明
	template <typename T>
		class GneVector;
	
	template <typename T>
		class ConstVectorView;
	
	template <typename T>
		class AlmostVector;
	
	
	template <typename T>
		class VectorView;
	
	template <typename T>
		class Vector;
	//
	// 向量赋值
	//
	template<typename T>
		class CommaAssignmentForVector
	{
		public:
		CommaAssignmentForVector( AlmostVector<T> &m, const T &t)
		{
#if NSL_RANGE_CHECK
			if ( m.size() == 0) 
			{
				throw gsl_error( "CommaAssignmentForVector: Vector has not been sized.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			m_    = &m;
			index = 0;
			m_->ref(index) = t;
		}
		
		CommaAssignmentForVector& operator, ( const T &t)
		{
			++index;
#if NSL_RANGE_CHECK
			if ( index >= m_->size()) 
			{
				throw gsl_error( "CommaAssignmentForVector::operator,(t): Assignment outside Vector range.",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			m_->ref(index) = t;
			return *this;
		}
		private:
		AlmostVector<T> *m_;
		size_t index;
	};
	
	//基类
	template <typename T>
		class GneVector 
	{
#ifndef NSL_SUPPORT_MULTITHREAD
	protected:
		
		static Vector<T> cacheVector0_;
		static Vector<T> cacheVector1_;
		static Vector<T> cacheVector2_;
#endif
	public:
		
		//constructor
		GneVector()
		{}
		virtual ~GneVector()
		{}
		
		// this converts the Vector , VectorView, and the ConstVectorView to valarray
		//operator valarray<T>()	const { return constVal() ;}
		
		//三个类都要重写的方法
		virtual const slice& getSlice() const = 0;
		virtual const valarray<T>& constVal() const = 0;
		
		//三个类公有的方法
		size_t size() 			const { return getSlice().size() ;}
		size_t start()			const {	return getSlice().start() ;}
		size_t stride() 		const { return getSlice().stride() ;}
		size_t index(size_t i)  const { return getSlice().start() + i*getSlice().stride() ;}
		
		size_t typesize()		const { return sizeof(T) ;}
		
		valarray<T> toVal()		const { return constVal()[getSlice()] ;} 
		
		
		inline T cref(size_t i) const 
		{ 
#if NSL_RANGE_CHECK
			if (i >= size())
			{
				gsl_error("Matrix::cref(i), index i out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return constVal()[index(i)] ;
		} 
		
		//Unary operators
		Vector<T> operator+ () const 
		{
			return *this;
		}
		
		Vector<T> operator- () const 
		{
			Vector<T> cache( *this, getSlice() );
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				cache[i] *= T(-1);
			}
			return cache;
		}
		
		T sum()		const 
		{ 
			T _sum = T(0);
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				_sum += cref(i);
			}
			
			return _sum;
		}
		
		T sumAbs() 	const 
		{ 
			T _sum = T(0);
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				_sum += nMath::abs(cref(i));
			}
			
			return _sum;
		}
		
		T min(size_t* index_out = NULL)		const 
		{ 
			size_t _index = 0;
			T _min = cref(0);
			size_t _size = size();
			for ( size_t i = 1; i < _size; ++i) 
			{
				if ( cref(i) < _min) 
				{
					_min  = cref(i);
					_index = i;
				}
			}
			
			if (index_out)
			{
				*index_out = _index;
			}
			
			return _min;
		}
		T max(size_t* index_out = NULL)		const 
		{ 
			size_t _index = 0;
			T _max = cref(0);
			size_t _size = size();
			for ( size_t i = 1; i < _size; ++i) 
			{
				if ( cref(i) > _max) 
				{
					_max  = cref(i);
					_index = i;
				}
			}
			
			if (index_out)
			{
				*index_out = _index;
			}
			
			return _max;
		}
		
		
		T minAbs(size_t* index_out = NULL)	const 
		{ 
			size_t _index = 0;
			T _min = cref(0);
			size_t _size = size();
			for ( size_t i = 1; i < _size; ++i) 
			{
				if ( abs(cref(i)) < abs(_min)) 
				{
					_min  = cref(i);
					_index = i;
				}
			}
			
			if (index_out)
			{
				*index_out = _index;
			}
			
			return _min;
		}
		
		T maxAbs(size_t* index_out = NULL)	const 
		{ 
			size_t _index = 0;
			T _max = cref(0);
			size_t _size = size();
			for ( size_t i = 1; i < _size; ++i) 
			{
				if ( abs(cref(i)) > abs(_max)) 
				{
					_max  = cref(i);
					_index = i;
				}
			}
			
			if (index_out)
			{
				*index_out = _index;
			}
			
			return _max;
		}
		
		
		T norm()	const 
		{ 
			T cache = 0;
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				cache += cref(i) * cref(i);
			}
			return ::sqrt(cache);
		}
		
		
		T normSq(const T scale = T(1)) 	const 
		{ 
			T cache = 0;
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				cache += (scale*cref(i)) * (scale*cref(i));
			}
			return ::sqrt(cache);
		}
		
		T average()	const 
		{
			return sum()/double(size());
		}
		
		Vector<T> unit()	const 
		{
			Vector<T> cache(*this, getSlice());
			T _norm = cache.norm();
			size_t _size = cache.size();
			for (size_t i=0; i<_size; i++)
			{
				cache[i] /= _norm;
			}
			return cache;
		}
		
		bool hasSubVector(const slice & s)	const
		{
			if (s.start() < start() 
				|| s.start() >= size()
				|| s.stride() >= size()
				|| s.size() > size())
			{
				return false;
			}
			return true;
		}
		
		bool isZero()	const 
		{ 
			size_t _size = size();
			for (size_t i=0; i<_size; i++)
			{
				if (nMath::abs(cref(i)) > nMath::ZERO)
				{
					return false;
				}
			}
			return true;
		}
		
		
		//
		// Binary Operators with type T.
		//
		friend Vector<T> operator* ( const GneVector<T> & gv, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) * t;
			}
			return cacheVector0_;
		}
		friend Vector<T> operator/ ( const GneVector<T> & gv, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) / t;
			}
			return cacheVector0_;
		}
		friend Vector<T> operator+ ( const GneVector<T> & gv, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) + t;
			}
			return cacheVector0_;;
		}
		friend Vector<T> operator- ( const GneVector<T> & gv, const T &t)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) - t;
			}
			return cacheVector0_;
		}
		friend Vector<T> operator* ( const T &t, const GneVector<T> & gv)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) * t;
			}
			return cacheVector0_;
		}
		friend Vector<T> operator/ ( const T &t, const GneVector<T> & gv)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = t / gv.cref(i);
			}
			return cacheVector0_;
		}
		friend Vector<T> operator+ ( const T &t, const GneVector<T> & gv)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = gv.cref(i) + t;
			}
			return cacheVector0_;
		}
		friend Vector<T> operator- ( const T &t, const GneVector<T> & gv)
		{
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(gv.size());
			for (size_t i=0; i<gv.size(); i++)
			{
				cacheVector0_(i) = t - gv.cref(i);
			}
			return cacheVector0_;
		}
		
		//
		// Binary Operators with type valarray<T>.
		//
		
		friend Vector<T> operator* ( const GneVector & v1, const GneVector<T> &v2)
		{
#if NSL_RANGE_CHECK
			if ( v1.size() != v2.size()) 
			{
				throw gsl_error( "valarray<T> operator*: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(v1.size());
			for (size_t i=0; i<v1.size(); i++)
			{
				cacheVector0_(i) = v1.cref(i) * v2.cref(i);
			}
			return cacheVector0_;
		}
		
		
		friend Vector<T> operator/ ( const GneVector & v1, const GneVector<T> &v2)
		{
#if NSL_RANGE_CHECK
			if ( v1.size() != v2.size()) 
			{
				throw gsl_error( "valarray<T> operator/: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(v1.size());
			for (size_t i=0; i<v1.size(); i++)
			{
				cacheVector0_(i) = v1.cref(i) / v2.cref(i);
			}
			return cacheVector0_;
		}
		
		
		friend Vector<T> operator+ ( const GneVector & v1, const GneVector<T> &v2)
		{
#if NSL_RANGE_CHECK
			if ( v1.size() != v2.size()) 
			{
				throw gsl_error( "valarray<T> operator+: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(v1.size());
			for (size_t i=0; i<v1.size(); i++)
			{
				cacheVector0_(i) = v1.cref(i) + v2.cref(i);
			}
			return cacheVector0_;
		}
		
		
		friend Vector<T> operator- ( const GneVector & v1, const GneVector<T> &v2)
		{
#if NSL_RANGE_CHECK
			if ( v1.size() != v2.size()) 
			{
				throw gsl_error( "valarray<T> operator-: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(v1.size());
			for (size_t i=0; i<v1.size(); i++)
			{
				cacheVector0_(i) = v1.cref(i) - v2.cref(i);
			}
			return cacheVector0_;
		}
		
		
		friend bool operator == (const GneVector<T>& _L, const GneVector<T>& _R)
		{
			size_t _size = _L.size();
			if (_size != _R.size())
			{
				return false;
			}
			for (size_t i=0; i<_size; i++)
			{
				if (_L.cref(i) != _R.cref(i))
				{
					return false;
				}
			}
			return true;
			
		}
		
		friend bool operator != (const GneVector<T>& _L, const GneVector<T>& _R)
		{
			return !(_L == _R);
			
		}
		
		//friend function
		
		friend void swap( AlmostVector<T> & avec1,  AlmostVector<T> & avec2)
		{
#if NSL_RANGE_CHECK
			if (avec1.size() != avec2.size())
			{
				throw gsl_error("swap function error:AlmostVector<T> size not compatible ",__FILE__,
					__LINE__, GSL_EBADLEN);
			}
#endif	
#ifdef NSL_SUPPORT_MULTITHREAD
			Vector<T> cacheVector0_;
#endif
			cacheVector0_.resize(avec1.size());
			assign(avec1,cacheVector0_);
			assign(avec2, avec1);
			assign(cacheVector0_, avec2);
		}
		
		friend std::ostream& operator<< ( std::ostream &os, const GneVector<T> &avec)
		{
			os << avec.constVal()[avec.getSlice()] << std::endl;
			return os;
		}
		
}; // end of class GneVector<T>

template <typename T>
class ConstVectorView : public GneVector<T> //常视图
{
	private:
		const valarray<T>* pVectorData_;
		slice viewSlice_;
		
		
	public:
		
		//重写基类的虚函数
		virtual const slice& getSlice() 			const { return viewSlice_ ;}
		virtual const valarray<T>& constVal() 		const { return *pVectorData_ ;}
		
		
		//constructor
		ConstVectorView()											:pVectorData_(NULL), viewSlice_(0,0,0)
		{}
		
		virtual ~ConstVectorView()									
		{ const_cast<const valarray<T>*>(pVectorData_) = NULL;}
		
		ConstVectorView(const valarray<T> & val)					:pVectorData_(&val), viewSlice_(0, val.size(), 1)
		{}
		
		ConstVectorView(const valarray<T> & val, const slice& s)	:pVectorData_(&val), viewSlice_(s)
		{}
		
		
		//以下两个构造函数是为了 
		//Vector			 ---> ConstVectorView
		//VectorView		 ---> ConstVectorView
		//ConstVectorView	 ---> ConstVectorView
		ConstVectorView(const GneVector<T> & vec)					:pVectorData_(&vec.constVal()), viewSlice_(vec.getSlice())
		{}
		
		ConstVectorView(const GneVector<T> & vec, const slice& s)	:pVectorData_(&vec.constVal()), viewSlice_(s)
		{}
		
		//subVector
		ConstVectorView<T> subVector(const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		ConstVectorView<T> operator [](const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		ConstVectorView<T> operator ()(const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		
		//access function
		T operator [](size_t i) const { return cref(i) ;}
		T operator ()(size_t i) const { return cref(i) ;}
		
}; // end of class ConstVectorView<T>


//VectorView and Vector should inherit this class
template <typename T>
class AlmostVector : public GneVector<T>
{
	
	public:
		AlmostVector()
		{}
		
		virtual ~AlmostVector()
		{}
		
		//using GneVector<T>::getSlice;
		
		//override this function in VectorView and Vector
		virtual valarray<T>& val() = 0;
		
		//VectorView and Vector 的共有函数
		
		//subVector
		ConstVectorView<T> subVector(const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		ConstVectorView<T> operator [](const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		ConstVectorView<T> operator ()(const slice & s)	const {return ConstVectorView<T>(*this, s) ;}
		
		VectorView<T> subVector(const slice & s)						{ return VectorView<T>(*this, s) ;}	
		VectorView<T> operator[] ( const slice & s)						{ return VectorView<T>(*this, s) ;}
		VectorView<T> operator() ( const slice & s)						{ return VectorView<T>(*this, s) ;}
		
		
		inline T& ref(size_t i)  
		{ 
#if NSL_RANGE_CHECK
			if (i >= size())
			{
				gsl_error("Matrix::ref(i), index i out of range",
					__FILE__, __LINE__, GSL_EBADLEN);
			}
#endif
			return val()[index(i)] ;
		} 
		
		
		
		//access function
		T& operator [](size_t i)  { return ref(i) ;}
		T& operator ()(size_t i)  { return ref(i) ;}
		
		T operator [](size_t i) const { return cref(i) ;}
		T operator ()(size_t i) const { return cref(i) ;}
		
		//make all elements zero
		void zero()				  { val() = T(0);}
		
		void setAll(const T & t)
		{
			val()[getSlice()] = T(t);
		}
		void swapElements(size_t i, size_t j)	
		{
			T cache = ref(i);
			ref(i) = ref(j);
			ref(j) = cache;
		}
		
		void setBasis(size_t i)
		{
			val()[getSlice()] = T(0);
			ref(i) = T(1);
		}
		
		void reverse()
		{
			size_t i=0;
			size_t _size = size();
			size_t m = _size/2;
			for (i=0; i<m; i++)
			{
				size_t j = _size - i - 1 ;
				T cache = ref(j);
				ref(j) = ref(i);
				ref(i) = cache;
			}
		}
		
		//compute and assign function
		void operator*= (const GneVector<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator*=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] *= v.constVal()[v.getSlice()];
		}
		
		void operator/= (const GneVector<T> &v) 
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator/=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] /= v.constVal()[v.getSlice()];
		}
		
		void operator+= (const GneVector<T> &v) 
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator+=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] += v.constVal()[v.getSlice()];
		}
		
		void operator-= (const GneVector<T> &v) 
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator-=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] -= v.constVal()[v.getSlice()];
		}
		
		//
		// Computed assignment operators.
		//
		void operator*= (const T &t)
		{
			for (size_t i=0; i<size(); i++)
			{
				ref(i) *= t;
			}
		}
		void operator/= (const T &t)
		{
			for (size_t i=0; i<size(); i++)
			{
				ref(i) /= t;
			}
		}
		void operator+= (const T &t)
		{
			for (size_t i=0; i<size(); i++)
			{
				ref(i) += t;
			}
		}
		void operator-= (const T &t)
		{
			for (size_t i=0; i<size(); i++)
			{
				ref(i) -= t;
			}
		}
		
		
		
		
		CommaAssignmentForVector<T> operator<< ( const T &t)
		{
			return CommaAssignmentForVector<T>( *this, t);
		}
		
		//
		// IO Stream functions.
		//
		
		friend std::istream& operator>> ( std::istream &is, AlmostVector<T> &avec)
		{
			for ( size_t i = 0; i < avec.size(); ++i) 
			{
				is >> avec[i];
			}
			return is;
		}
};

template <typename T>
class VectorView : public AlmostVector<T> //可变视图
{
	private:
		valarray<T>* pVectorData_;
		slice viewSlice_;
		
		
	public:
		
		//重写基类的虚函数
		virtual const slice& getSlice() 			const { return viewSlice_ ;}
		virtual const valarray<T>& constVal() 		const { return *pVectorData_ ;}
		
		//constructor
		VectorView()												:pVectorData_(NULL), viewSlice_(0,0,0)
		{}
		
		virtual ~VectorView()								
		{ const_cast<valarray<T>*>(pVectorData_) = NULL ;}
		
		//下面两个构造函数是为了
		//Vector			 ---> MatrixView
		//VectorView		 ---> VectorView
		VectorView(AlmostVector<T> & vec)							:pVectorData_(&vec.val()), viewSlice_(vec.getSlice())
		{}
		
		VectorView(AlmostVector<T> & vec, const slice& s)			:pVectorData_(&vec.val()), viewSlice_(s)
		{}
		
		VectorView(valarray<T> & val)						:pVectorData_(&val), viewSlice_(0, val.size(), 1)
		{
		}
		
		VectorView(valarray<T> & val, const slice& s)		:pVectorData_(&val), viewSlice_(s)
		{}
		
		
		virtual valarray<T>& val() 	 { return *pVectorData_ ;}
		
		VectorView<T> operator= ( const AlmostVector<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] = v.constVal()[v.getSlice()];
			
			return *this;
		}
		
		VectorView<T> operator= ( const valarray<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator=: valarray not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] = v;
			return *this;
		}
		
		// This is for A[i] = vector.
		VectorView<T> operator= ( const std::vector<T> &v)
		{
#if NSL_RANGE_CHECK
			if ( size() != v.size()) 
			{
				throw gsl_error( "void operator=: vector not compatible.", __FILE__, __LINE__, GSL_EINVAL);
			}
#endif
			val()[getSlice()] = getValArrayFromVector( v);
			return *this;
		}
		
		VectorView<T> operator= ( const T & t)
		{
			val()[getSlice()] = t;
			return *this;
		}
		
}; // end of class VectorView<T>




//向量类
template <typename T>
class Vector : public AlmostVector<T>
{
	private:
		valarray<T> vectorData_;
		slice		slice_;
		
	public:
		
		//重写基类的虚函数
		virtual const slice& getSlice() 			const   { return slice_ ;}
		virtual const valarray<T>& constVal() 		const   { return vectorData_ ;}
		
		//another way to get the valarray<T>&
		virtual valarray<T>& val()									{ return vectorData_ ;}
		
		//destruction
		virtual ~Vector()
		{}
		
		//cosntruction
		Vector()											: vectorData_(0), slice_(0, 0, 0)
		{}
		
		//由STL valarray 构造
		Vector(const valarray<T>& val) 						: vectorData_(val), slice_(0, val.size(), 1)
		{}
		
		Vector(const valarray<T>& val, const slice& s) 		: vectorData_(val), slice_(s)
		{}
		
		
		//给定大小
		Vector(size_t n) 									: vectorData_(T(0), n), slice_(0, n, 1)
		{}
		
		//给定大小和值
		Vector(size_t n, const T & val) 					: vectorData_(val, n), slice_(0, n, 1)
		{}
		
		
		//由数组构造
		Vector(const T * p, const slice & s) 				: vectorData_(p, s.size()), slice_(s)
		{}
		
		Vector(const T * p, size_t n) 						: vectorData_(p, n), slice_(0, n, 1)
		{}
		
		
		//以下两个构造函数是为了 
		//Vector			 ---> Vector
		//VectorView		 ---> Vector
		//ConstVectorView	 ---> Vector
		Vector(const GneVector<T> & vec) 				: vectorData_( (vec.constVal()[vec.getSlice()])), slice_(0, vec.size(), 1)
		{}
		
		Vector(const GneVector<T> & vec, const slice & s) : vectorData_( (vec.constVal())[s]), slice_(s)
		{}
		
		
		void resize(size_t n)	{ vectorData_.resize(n); slice_ = slice(0, n, 1) ;}
		
		
		Vector<T>& operator = (const GneVector<T> & other)
		{
			if (size() != other.size())
			{
				resize(other.size());
			}
			vectorData_ = other.constVal()[other.getSlice()];
			slice_ = slice(0, other.size(), 1);
			
			return *this;
		}
		
		Vector<T>& operator= ( const valarray<T> &v)
		{
			val() = v;
			slice_ = slice(0, v.size(), 1);
			return *this;
		}
		
		// This is for A[i] = vector.
		Vector<T>& operator= ( const std::vector<T> &v)
		{
			val() = getValArrayFromVector( v);
			slice_ = slice(0, v.size(), 1);
			return *this;
		}
		
		Vector<T> operator= ( const T & t)
		{
			val() = t;
			slice_ = slice(0, size(), 1);
			return *this;
		}
};  // end of class Vector<T>

#ifndef NSL_SUPPORT_MULTITHREAD
// Vector evaluation caches.
template <typename T> Vector<T> GneVector<T>::cacheVector0_;
template <typename T> Vector<T> GneVector<T>::cacheVector1_;
template <typename T> Vector<T> GneVector<T>::cacheVector2_;
#endif

#ifndef __valarray_functions_H
#define __valarray_functions_H
//
// Functions for valarray.
//
template <typename T> inline void assign(const AlmostVector<T>& v1, AlmostVector<T>& v2)
{
	for (size_t i=0; i<v1.size(); i++)
	{
		v2(i) = v1(i);
	}
}

template <typename T> inline size_t absMin( const GneVector<T> &v, size_t start = 0)
{
#if NSL_RANGE_CHECK
    if ( v.size() == 0) 
	{
		throw gsl_error( "size_t absMin(): Vector is empty.", __FILE__, __LINE__, GSL_EINVAL);
    }
#endif
    size_t index = start;
    T dMin = abs( v.cref(index));
    for ( size_t i = start + 1; i < v.size(); ++i) 
	{
		if ( abs( v.cref(i)) < dMin) 
		{
			dMin  = v.cref(i);
			index = i;
		}
    }
    return index;
}

template <typename T> inline size_t absMax( const GneVector<T> &v, size_t start = 0)
{
#if NSL_RANGE_CHECK
	if ( v.size() == 0) 
	{
		throw gsl_error( "size_t absMax(): Vector is empty.", __FILE__, __LINE__, GSL_EINVAL);
	}
#endif
    size_t index = start;
    T dMin = abs( v.cref(index));
    for ( size_t i = start + 1; i < v.size(); ++i) 
	{
		if ( abs( v.cref(i)) > dMin) 
		{
			dMin  = v.cref(i);
			index = i;
		}
    }
    return index;
}

template <typename T> inline size_t min( const GneVector<T> &v, size_t start = 0)
{
#if NSL_RANGE_CHECK
	if ( v.size() == 0) 
	{
		throw gsl_error( "size_t min(): Vector is empty.", __FILE__, __LINE__, GSL_EINVAL);
	}
#endif
	size_t index = start;
	T dMin = v.cref(index);
	for ( size_t i = start + 1; i < v.size(); ++i) 
	{
		if ( v.cref(i) < dMin) {
			dMin  = v.cref(i);
			index = i;
		}
	}
    return index;
}

template <typename T> inline size_t max( const GneVector<T> &v, size_t start = 0)
{
#if NSL_RANGE_CHECK
    if ( v.size() == 0) 
	{
		throw gsl_error( "size_t max(): Vector is empty.", __FILE__, __LINE__, GSL_EINVAL);
    }
#endif
    size_t index = start;
    T dMin = v.cref(index);
    for ( size_t i = start + 1; i < v.size(); ++i) 
	{
		if ( v.cref(i) > dMin) 
		{
			dMin  = v.cref(i);
			index = i;
		}
    }
    return index;
}

//
// The following dot product functions are duplicated for speed increase.
//
template <typename T> inline T dot( const GneVector<T> &v1, const GneVector<T> &v2)
{
#if NSL_RANGE_CHECK
    if ( v1.size() != v2.size()) 
	{
		throw gsl_error( "dot(v1,v2): Vectors must be same length.", __FILE__, __LINE__, GSL_EINVAL);
    }
#endif
    T dp = T(0);
    for ( size_t i = 0; i < v1.size(); i++)
	{
		dp += v1.cref(i) * v2.cref(i);
    }
    return dp;
}

//
// The following norm functions are duplicated for speed increase.
//
template <typename T> inline T norm( const GneVector<T> &v1)
{
    T t(0);
    for ( size_t i = 0; i < v1.size(); ++i) 
	{
		t += v1.cref(i) * v1.cref(i);
    }
    return ::sqrt( t);
}

template <typename T> inline bool isZero( const GneVector<T> &v1)
{
    for ( size_t i = 0; i < v1.size(); ++i) 
	{
		if ( nMath::abs( v1.cref(i)) > nMath::ZERO) 
		{
			return false;
		}
    }
    return true;
}

//
// IO Stream functions for valarray
//
template <typename T> std::ostream& operator<< ( std::ostream &os, const valarray<T> &v)
{
	if ( v.size() == 0) 
	{
		os << "Null Vector";
	}
	else 
	{
		os.width(11);
		if ( nMath::abs( v[0]) < nMath::ZERO) 
		{
			os << 0;
		}
		else {
			os << v[0];
		}
		for ( size_t i=1; i < v.size(); ++i) 
		{
			os << ' ';
			os.width(11);
			if ( nMath::abs( v[i]) < nMath::ZERO) 
			{
				os << 0;
			}
			else 
			{
				os << v[i];
			}
		}
	}
	return os;
}

template <typename T> std::istream& operator>> ( std::istream &is, valarray<T> &v)
{
	for ( size_t i = 0; i < v.size(); ++i) 
	{
		is >> v[i];
	}
	return is;
}
#endif // __valarray_functions_H


} // end of namespace gslcpp

#endif // NSL_VECTOR_H__
