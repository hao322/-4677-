#pragma once
#include<gdal_priv.h>
#include<iostream>

/*
分别定义三个数组存储R、G、B
再分别定义三个数组存储H、S、I
中间一个临时数组存储theta
*/



class HIS_ImgFusion
{
public:
	HIS_ImgFusion(const char*);//有参构造，传入影像的地址进行初始化
	~HIS_ImgFusion();
	int getXsize();
	int getBandnum();
	int getYsize();

private:
	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;

};



void Tttt(const char*, const char*);