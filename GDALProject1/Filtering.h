#pragma once
#include <iostream>
#include <gdal_priv.h>


class Filtering
{
public:
	Filtering(const char*);//构造函数
	int getXsize();
	int getBandnum();
	unsigned char* Mean();
	unsigned char* Median(); //返回处理后的图像数据
	void SaveFile(const char*, unsigned char*);
private:
	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;
};

