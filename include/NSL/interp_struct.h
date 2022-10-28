/********************************************************************
	filename: 	interp_struct.h
	author:		hu zhijian
	created:	5:5:2010   16:57
	brief:	提供给所有插值类的接口数据类型
*********************************************************************/

#ifndef NWPU_INTRP_STRUCT_H
#define NWPU_INTRP_STRUCT_H

typedef struct
{
	size_t	 dimension; //插值维数

	double** x_Array;	//自变量

	size_t*  x_Array_Num; //每一维自变量的个数

	double*	 y_Array;	//应变量

	size_t   y_Array_Num; //应变量的个数

	bool	 extro_interp_flag; //=true时可外插, 否则出错提示

	size_t	 Lagrange_point_Num; //拉格朗日插值点数
	//在拉格朗日插值中 point_Num =2时为线性插值，
	//=3时为抛物线插值 
	//point_Num越大，插值结果越精确，但是速度也越慢,有可能发生龙格现象
						//建议小于5，大于1
}interp_table;


#endif