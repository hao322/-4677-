#pragma once
#include <iostream>
#include <gdal_priv.h>

class Data
{
public:
	Data(const char*);//�вι��죬����Ӱ��ĵ�ַ���г�ʼ��
	int getXsize();
	int getBandnum();
	int getYsize();
	std::pair<int, int> match_ncc(Data &);
	float cal_Mean(unsigned char*, int);
	float cal_NCC(unsigned char*, unsigned char*, int);

private:
	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;
};

