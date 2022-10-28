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

#ifndef NSL_VAR

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


/* Turn range checking on by default, unless the user defines
NSL_RANGE_CHECK_OFF, or defines GSL_RANGE_CHECK to 0 explicitly */

#ifdef NSL_RANGE_CHECK_OFF
# ifndef NSL_RANGE_CHECK
#  define NSL_RANGE_CHECK 0
# else
#  error "cannot set both GSL_RANGE_CHECK and GSL_RANGE_CHECK_OFF"
# endif
#else
# ifndef NSL_RANGE_CHECK
#  define NSL_RANGE_CHECK 1
# endif
#endif


#define NSL_SUPPORT_MULTITHREAD

#endif // NSL_NSL_H__