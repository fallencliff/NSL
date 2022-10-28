# Microsoft Developer Studio Project File - Name="NSL" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Dynamic-Link Library" 0x0102
# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=NSL - Win32 Debug
!MESSAGE This is not a valid makefile. To build this project using NMAKE,
!MESSAGE use the Export Makefile command and run
!MESSAGE 
!MESSAGE NMAKE /f "NSL.mak".
!MESSAGE 
!MESSAGE You can specify a configuration when running NMAKE
!MESSAGE by defining the macro CFG on the command line. For example:
!MESSAGE 
!MESSAGE NMAKE /f "NSL.mak" CFG="NSL - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "NSL - Win32 Release" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "NSL - Win32 Debug" (based on "Win32 (x86) Dynamic-Link Library")
!MESSAGE "NSL - Win32 StaticDebug" (based on "Win32 (x86) Static Library")
!MESSAGE "NSL - Win32 StaticRelease" (based on "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""

!IF  "$(CFG)" == "NSL - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\build\win32\dynamic\Release"
# PROP Intermediate_Dir "obj"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
CPP=cl.exe
# ADD BASE CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "NSL_EXPORTS" /Yu"stdafx.h" /FD /c
# ADD CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "NSL_DLL" /D "DLL_EXPORT" /FR /FD /c
# SUBTRACT CPP /YX /Yc /Yu
MTL=midl.exe
# ADD BASE MTL /nologo /D "NDEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "NDEBUG" /mktyplib203 /win32
RSC=rc.exe
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo /o"obj/NSL.bsc"
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /machine:I386
# ADD LINK32 /nologo /dll /pdb:"obj/NSL.pdb" /machine:I386
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\build\win32\dynamic\Debug"
# PROP Intermediate_Dir "obj"
# PROP Ignore_Export_Lib 0
# PROP Target_Dir ""
CPP=cl.exe
# ADD BASE CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "NSL_EXPORTS" /Yu"stdafx.h" /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "_MBCS" /D "_USRDLL" /D "NSL_DLL" /D "DLL_EXPORT" /FD /GZ /c
# SUBTRACT CPP /YX /Yc /Yu
MTL=midl.exe
# ADD BASE MTL /nologo /D "_DEBUG" /mktyplib203 /win32
# ADD MTL /nologo /D "_DEBUG" /mktyplib203 /win32
RSC=rc.exe
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib user32.lib gdi32.lib winspool.lib comdlg32.lib advapi32.lib shell32.lib ole32.lib oleaut32.lib uuid.lib odbc32.lib odbccp32.lib /nologo /dll /debug /machine:I386 /pdbtype:sept
# ADD LINK32 /nologo /dll /pdb:"obj/NSL.pdb" /debug /machine:I386 /pdbtype:sept
# SUBTRACT LINK32 /pdb:none

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "..\build\win32\static\debug"
# PROP Intermediate_Dir "obj"
# PROP Target_Dir ""
MTL=midl.exe
CPP=cl.exe
# ADD BASE CPP /nologo /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /YX /FD /GZ /c
# ADD CPP /nologo /MTd /W3 /Gm /GX /ZI /Od /D "WIN32" /D "_DEBUG" /D "_MBCS" /D "_LIB" /D "NSL_DLL" /D "DLL_EXPORT" /YX /FD /GZ /c
RSC=rc.exe
# ADD BASE RSC /l 0x804 /d "_DEBUG"
# ADD RSC /l 0x804 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "..\build\win32\static\release"
# PROP Intermediate_Dir "obj"
# PROP Target_Dir ""
MTL=midl.exe
CPP=cl.exe
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /YX /FD /c
# ADD CPP /nologo /MT /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_MBCS" /D "_LIB" /D "NSL_DLL" /D "DLL_EXPORT" /YX /FD /c
RSC=rc.exe
# ADD BASE RSC /l 0x804 /d "NDEBUG"
# ADD RSC /l 0x804 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo

!ENDIF 

# Begin Target

# Name "NSL - Win32 Release"
# Name "NSL - Win32 Debug"
# Name "NSL - Win32 StaticDebug"
# Name "NSL - Win32 StaticRelease"
# Begin Group "Source Files"

# PROP Default_Filter "cpp;c;cxx;rc;def;r;odl;idl;hpj;bat"
# Begin Group "block"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\block\Block.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# End Group
# Begin Group "combination"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\combination\Combination.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# End Group
# Begin Group "complex"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\complex\Complex.cpp
# End Source File
# End Group
# Begin Group "interpolation"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\interpolation\baseClass\BaseInterp.cpp
# End Source File
# Begin Source File

SOURCE=..\src\interpolation\BiCubicSpline\BiCubicSpline.cpp
# End Source File
# Begin Source File

SOURCE=..\src\interpolation\cubicSpline\CubicSpline.cpp
# End Source File
# Begin Source File

SOURCE=..\src\interpolation\PolyCubicSpline\PolyCubicSpline.cpp
# End Source File
# Begin Source File

SOURCE=..\src\interpolation\PolyLagrangeInterp\PolyLagrangeInterp.cpp
# End Source File
# End Group
# Begin Group "matrix"

# PROP Default_Filter "cpp"
# End Group
# Begin Group "vector"

# PROP Default_Filter "cpp;h"
# Begin Source File

SOURCE=..\include\NSL\Vector.h
# End Source File
# End Group
# Begin Group "random"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\random\RandGenerator.cpp
# End Source File
# End Group
# Begin Group "randist"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\randist\Randist.cpp
# End Source File
# End Group
# Begin Group "permutation"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\permutation\Permutation.cpp
# End Source File
# End Group
# Begin Group "montecarlo"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\MonteCarlo\MonteCarlo.cpp
# End Source File
# End Group
# Begin Group "cmath"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\src\mathconst\MathConst.cpp
# End Source File
# End Group
# Begin Group "error"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\err\error.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\err\error_message.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\err\error_stream.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=..\src\err\error_strerror.cpp

!IF  "$(CFG)" == "NSL - Win32 Release"

# SUBTRACT CPP /YX /Yc /Yu

!ELSEIF  "$(CFG)" == "NSL - Win32 Debug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticDebug"

!ELSEIF  "$(CFG)" == "NSL - Win32 StaticRelease"

!ENDIF 

# End Source File
# End Group
# Begin Group "histogram"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\histogram\histogram.cpp
# End Source File
# Begin Source File

SOURCE=..\src\histogram\histogram2D.cpp
# End Source File
# End Group
# Begin Group "statistics"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE=..\src\statistics\Statis.cpp
# End Source File
# End Group
# Begin Group "ode-initval"

# PROP Default_Filter "cpp"
# Begin Source File

SOURCE="..\src\ode-initval\Odeiv.cpp"
# End Source File
# End Group
# Begin Group "qrng"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\src\qrng\Qrng.cpp
# End Source File
# End Group
# Begin Group "ntuple"

# PROP Default_Filter ""
# Begin Source File

SOURCE=..\src\ntuple\Ntuple.cpp
# End Source File
# End Group
# Begin Group "cblas"

# PROP Default_Filter "cpp"
# End Group
# End Group
# Begin Group "Header Files"

# PROP Default_Filter "h;hpp;hxx;hm;inl"
# Begin Source File

SOURCE=..\include\NSL\Block.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Combination.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Complex.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\gsl_errno.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\histogram.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\MathConst.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Matrix.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\MonteCarlo.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\NSL.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Ntuple.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Odeiv.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Permutation.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\PolyCubicSpline.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\PolyLagrangeInterp.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Qrng.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\RandGenerator.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Randist.h
# End Source File
# Begin Source File

SOURCE=..\include\NSL\Statis.h
# End Source File
# End Group
# End Target
# End Project
