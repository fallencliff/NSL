/**
 *@file GBlock.h
 *@brief
 *@author
 *@date 2010-03-21
 */

#ifndef NSL_BLOCK_H__
#define NSL_BLOCK_H__

#include <NSL.h>
#include <stdio.h>
#include <gsl_errno.h>
//#include <Vector.h>
#include <CommonStruct.h>
namespace gslcpp
{
	/**
	 * @brief CGBlock类
	 *
	 */
	
	class NSL_EXPORT CBlock
	{
		//friend class CVector;
		
		//friend class CMatrix;

		public:

		private:
			/** 
			 * @brief: 默认构造函数, 设为私有函数，不许使用
			 *
			 */
			CBlock();

		protected:

			gsl_block* gsldata;

			/** 
			 * @brief: 分配一块能存放n个双精度类型的内存, 内存不初始化,
			 *
			 */
			gsl_block* alloc(const size_t n);
			
			/** 
			 * @brief: 分配一块能存放n个双精度类型的内存, 内存初始化为0,
			 *
			 */
			gsl_block* calloc(const size_t n);

		public:
		
			/** 
			 * @brief: 带参构造函数, 
			 *
			 */
			CBlock(const size_t n, bool clear = true);
			
			/** 
			 * @brief: 拷贝构造函数, 
			 *
			 */
			CBlock(const CBlock& other);

			
			/** 
			 * @brief: 析构函数, 
			 *
			 */
			~CBlock();

			/** 
			 * @brief: 返回大小, 
			 *
			 * @org: 
			 */
			size_t Size() const;

			double * GetDataPtr() const;
			
			int Save(const char* filename, const char* format) const;

			int RawSave(const char* filename, 
				const size_t n,
				const size_t stride, 
                const char *format);
			int Load(const char* filename);
		
	};
}

#endif // NSL_BLOCK_H__