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
	 * @brief CGBlock��
	 *
	 */
	
	class NSL_EXPORT CBlock
	{
		//friend class CVector;
		
		//friend class CMatrix;

		public:

		private:
			/** 
			 * @brief: Ĭ�Ϲ��캯��, ��Ϊ˽�к���������ʹ��
			 *
			 */
			CBlock();

		protected:

			gsl_block* gsldata;

			/** 
			 * @brief: ����һ���ܴ��n��˫�������͵��ڴ�, �ڴ治��ʼ��,
			 *
			 */
			gsl_block* alloc(const size_t n);
			
			/** 
			 * @brief: ����һ���ܴ��n��˫�������͵��ڴ�, �ڴ��ʼ��Ϊ0,
			 *
			 */
			gsl_block* calloc(const size_t n);

		public:
		
			/** 
			 * @brief: ���ι��캯��, 
			 *
			 */
			CBlock(const size_t n, bool clear = true);
			
			/** 
			 * @brief: �������캯��, 
			 *
			 */
			CBlock(const CBlock& other);

			
			/** 
			 * @brief: ��������, 
			 *
			 */
			~CBlock();

			/** 
			 * @brief: ���ش�С, 
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