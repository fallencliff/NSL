/********************************************************************
	filename: 	blash_type.h
	author:		hu zhijian
	created:	16:5:2010   12:01
	brief:	
*********************************************************************/

#ifndef NSL_BLASH_TYPE_H__
#define NSL_BLASH_TYPE_H__




/*
 * Enumerated and derived types
 */
#define CBLAS_INDEX size_t  /* this may vary between platforms */

enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE {CblasLeft=141, CblasRight=142};


typedef  CBLAS_INDEX  CBLAS_INDEX_t;
typedef  enum CBLAS_ORDER       CBLAS_ORDER_t;
typedef  enum CBLAS_TRANSPOSE   CBLAS_TRANSPOSE_t;
typedef  enum CBLAS_UPLO        CBLAS_UPLO_t;
typedef  enum CBLAS_DIAG        CBLAS_DIAG_t;
typedef  enum CBLAS_SIDE        CBLAS_SIDE_t;


#endif // NSL_BLASH_TYPE_H__