/**
 *@file GBlock.h
 *@brief
 *@author
 *@date 2010-03-21
 */

#ifndef	GSLCPP_CGBLOCK_H
#define	GSLCPP_CGBLOCK_H

#include <stdio.h>
#include <gsl_errno.h>
#include <Vector.h>

#ifdef SL_DLL_EXPORTS
#define SL_DLL_API __declspec(dllexport)
#else
#define SL_DLL_API __declspec(dllimport)
#endif

namespace gslcpp
{
	/**
	 * @brief CGBlock��
	 *
	 */
	
	

	class SL_DLL_API CBlock
	{
		friend class CVector;
		
		friend class CMatrix;

		public:
			typedef struct
			{
				size_t size;
				double *data;
			
			}gsl_block;

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

			
			int Save(const char* filename, const char* format) const;

			int RawSave(const char* filename, 
				const size_t n,
				const size_t stride, 
                const char *format);
			int Load(const char* filename);
		
	};
}
#endif