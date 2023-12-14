#pragma once
#include <iostream>
#include <gdal_priv.h>

class Data
{
public:
	Data(const char*);//有参构造，传入影像的地址进行初始化
	int getXsize();
	int getBandnum();
	float* getNDVI();
	void saveFile(const char*, float*);

private:
	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;
};

