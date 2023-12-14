#pragma once
#include <gdal_priv.h>
#include <iostream>

class Initialize
{
public:
	Initialize(const char*); //有参构造
	~Initialize(); //析构，可设计一个保存输出的函数
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

//CVA变化检测
unsigned int* getChange(Initialize&, Initialize&);
//确定阈值
int threshold(unsigned int*, int, int);
//二值化
unsigned short* segmentation(unsigned int* change, int threshold, int Xsize, int Ysize);
//保存图像
void saveFile(unsigned short*, const char*, Initialize);

//精度评价
void EvaAccuracy(Initialize&, Initialize&);