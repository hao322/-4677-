#pragma once
#include <iostream>
#include <gdal_priv.h>

class Data
{
public:
	Data(const char*);//�вι��죬����Ӱ��ĵ�ַ���г�ʼ��
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

