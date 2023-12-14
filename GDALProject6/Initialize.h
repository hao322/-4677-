#pragma once
#include <gdal_priv.h>
#include <iostream>

class Initialize
{
public:
	Initialize(const char*); //�вι���
	~Initialize(); //�����������һ����������ĺ���
	int getXsize();
	int getYsize();
	int getBandnum();
	GDALDataset* getDataset();
	GDALDataType getDatatype();

private:
	GDALDataset* Dataset;
	GDALDataType Datatype;
	int Xsize;
	int Ysize;
	int Bandnum;

};

//CVA�仯���
unsigned int* getChange(Initialize&, Initialize&);
//ȷ����ֵ
int threshold(unsigned int*, int, int);
//��ֵ��
unsigned short* segmentation(unsigned int* change, int threshold, int Xsize, int Ysize);
//����ͼ��
void saveFile(unsigned short*, const char*, Initialize);

//��������
void EvaAccuracy(Initialize&, Initialize&);