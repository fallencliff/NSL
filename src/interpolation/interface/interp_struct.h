/********************************************************************
	filename: 	interp_struct.h
	author:		hu zhijian
	created:	5:5:2010   16:57
	brief:	�ṩ�����в�ֵ��Ľӿ���������
*********************************************************************/

#ifndef NWPU_INTRP_STRUCT_H
#define NWPU_INTRP_STRUCT_H

typedef struct
{
	size_t	 dimension; //��ֵά��

	double** x_Array;	//�Ա���

	size_t*  x_Array_Num; //ÿһά�Ա����ĸ���

	double*	 y_Array;	//Ӧ����

	size_t   y_Array_Num; //Ӧ�����ĸ���

	bool	 extro_interp_flag; //=trueʱ�����, ���������ʾ

	size_t	 Lagrange_point_Num; //�������ղ�ֵ����
	//���������ղ�ֵ�� point_Num =2ʱΪ���Բ�ֵ��
	//=3ʱΪ�����߲�ֵ 
	//point_NumԽ�󣬲�ֵ���Խ��ȷ�������ٶ�ҲԽ��,�п��ܷ�����������
						//����С��5������1
}interp_table;


#endif