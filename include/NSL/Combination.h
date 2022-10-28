/********************************************************************
	filename: 	Combination.h
	author:		hu zhijian
	created:	5:5:2010   16:50
	brief:	组合类
*********************************************************************/


#ifndef NSL_COMBINATION_H__
#define NSL_COMBINATION_H__


#include <NSL.h>
#include <CommonStruct.h>
namespace gslcpp
{
	class NSL_EXPORT CCombination
	{

	private:
		gsl_combination* gsldata;


	private:
		void Free();
		gsl_combination* Alloc(const size_t n, const size_t k);
		gsl_combination* Calloc(const size_t n, const size_t k);


	public:
		CCombination();
		virtual ~CCombination();
		
		CCombination(const size_t n, const size_t k, bool flag = false);
	
		void Resize(size_t n, size_t k, bool flag = false);
		
		bool IsValid();

		size_t* GetDataPtr()
		{
			return gsldata->data;
		}

		const size_t* GetDataPtr() const
		{
			return gsldata->data;
		}

		size_t Size()
		{
			return gsldata->k;
		}

		size_t GetRange()
		{
			return gsldata->n;
		}

		void InitFirst();
		void InitLast();
		
		bool Next();
		bool Prev();

		CCombination& operator=(const CCombination& other);
		size_t operator[](const size_t index) const;

		bool Save(const char* filename, const char* format) const;
		bool Load(const char* filename);

		void print(); //调试函数
	};

}


#endif // NSL_COMBINATION_H__