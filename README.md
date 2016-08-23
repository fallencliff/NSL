# GSL
The NSL is a independent c++ scientific library of GSL

declarations:
I didn't do anything except rewrite GSL in c++. so any version things i didn't mention applies with the GSL


使用NSL库

1.	把NSL include 头文件的文件夹设为VC的include path.
	
2.	链接选项中加入 NSL.lib(有静态库和动态库两种)

	若要使用NSL动态库
	在预处理器定义中加入NSL_DLL， 并把NSL.dll拷贝到应用程序可找到的路径内





if your want to rebuild the NSL in win32, follow the following steps
 
1.	set the include folder in your search path.

2.	set the src folder in your search path.

3.	use the vc6.0 project files in the project folder
