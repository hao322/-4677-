#pragma once
#include<gdal_priv.h>
#include<iostream>

/*
�ֱ�����������洢R��G��B
�ٷֱ�����������洢H��S��I
�м�һ����ʱ����洢theta
*/



class HIS_ImgFusion
{
public:
	HIS_ImgFusion(const char*);//�вι��죬����Ӱ��ĵ�ַ���г�ʼ��
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