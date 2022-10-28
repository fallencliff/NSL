/********************************************************************
	filename: 	NSL.h
	author:		hu zhijian
	created:	12:5:2010   17:12
	brief:	preprocessor definitions for export
*********************************************************************/

#ifndef NSL_NSL_H__
#define NSL_NSL_H__


// Insert your headers here
#define WIN32_LEAN_AND_MEAN		// Exclude rarely-used stuff from Windows headers

#ifndef GSL_VAR

#ifdef WIN32
#  ifdef NSL_DLL
#    ifdef DLL_EXPORT
#      define NSL_VAR extern __declspec(dllexport)
#      define NSL_EXPORT __declspec(dllexport)
#    else
#      define NSL_VAR extern __declspec(dllimport)
#      define NSL_EXPORT __declspec(dllimport)
#    endif
#  else
#    define NSL_VAR extern
#    define NSL_EXPORT
#  endif
#else
#  define NSL_VAR extern
#  define NSL_EXPORT
#endif

#endif


#endif // NSL_NSL_H__