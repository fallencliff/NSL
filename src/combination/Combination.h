#ifndef	NWPU_COMBINATION_H
#define NWPU_COMBINATION_H
/** 
 * @file Combination.h
 * @brief 组合类
 * @author 
 * @date 2010-04-13
 */

#include <SL_dll.h>

namespace gslcpp
{
	class SL_DLL_API CCombination
	{
	public:
		typedef struct
		{
			size_t n;
			size_t k;
			size_t *data;
		}gsl_combination;

	private:
		gsl_combination* gsldata;


	private:
		void Free();
		gsl_combination* Alloc(const size_t n, const size_t k);
		gsl_combination* Calloc(const size_t n, const size_t k);


	public:
		CCombination();
		~CCombination();
		
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


#endif