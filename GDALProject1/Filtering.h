#pragma once
#include <iostream>
#include <gdal_priv.h>


class Filtering
{
public:
	Filtering(const char*);//���캯��
	int getXsize();
	int getBandnum();
	unsigned char* Mean();
	unsigned char* Median(); //���ش�����ͼ������
	void SaveFile(const char*, unsigned char*);
private:
	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;
};

