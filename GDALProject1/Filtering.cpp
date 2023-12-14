#include "pch.h"
#include "Filtering.h"
#include <gdal_priv.h>
#include <iostream>
#include <algorithm>

using namespace std;

//有参构造函数赋初值
Filtering::Filtering(const char* datapath) {
	this->Dataset = (GDALDataset*)GDALOpen(datapath, GA_ReadOnly);
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
	this->DataType = Dataset->GetRasterBand(1)->GetRasterDataType();
}

int Filtering::getXsize()
{
	return this->Xsize;
}
int Filtering::getBandnum()
{
	return this->Bandnum;
}

//均值滤波
unsigned char* Filtering::Mean()
{
	//申请内存 一个波段
	unsigned char* Input = new unsigned char[Xsize * Ysize * DataType];
	unsigned char* Output = new unsigned char[Xsize * Ysize * DataType];

	GDALRasterBand* band1 = Dataset->GetRasterBand(1);
	band1->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Input, Xsize, Ysize, DataType, 0, 0);

	for (int i = 0; i < Ysize; i++) {
		for (int j = 0; j < Xsize; j++) {
			if (i == 0 || j == 0 || i == Ysize - 1 || j == Xsize - 1) {
				Output[Xsize * i + j] = Input[Xsize * i + j];
			}
		}
	}

	int temp[9] = { 0 };//3*3的
	for (int i = 1; i < Ysize - 1; i++) {
		for (int j = 1; j < Xsize - 1; j++) {
			for (int k = 0; k < 9; k++) {
				*(temp + k) = Input[Xsize * (i + k / 3) + j + k % 3];
			}
			int sum = 0;
			for (int t = 0; t < 9; t++) {
				sum += *(temp + t);
			}

			Output[Xsize * i + j] = sum / 9;
		}
	}
	return Output;
}

//中值滤波
unsigned char* Filtering::Median()
{
	unsigned char* Input = new unsigned char[Xsize * Ysize * DataType];
	unsigned char* Output = new unsigned char[Xsize * Ysize * DataType];

	GDALRasterBand* band1 = Dataset->GetRasterBand(1);
	band1->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Output, Xsize, Ysize, DataType, 0, 0);
	for (int i = 0; i < this->Ysize; i++) {
		for (int j = 0; j < this->Xsize; j++) {
			if (i == 0 || j == 0 || i == Ysize - 1 || j == Xsize - 1) {
				Output[Xsize * i + j] = Input[Xsize * i + j];
			}
		}
	}
	int temp[9] = { 0 };//3*3的
	for (int i = 1; i < this->Ysize - 1; i++) {
		for (int j = 1; j < this->Xsize - 1; j++) {
			for (int k = 0; k < 9; k++) {
				*(temp + k) = Input[Xsize * (i + (k / 3 - 1)) + j + k % 3 - 1];
			}
			sort(temp, temp + 9);
			Output[Xsize * i + j] = (unsigned char)*(temp + 4);
		}
	}

	return Output;
}

//保存影像
void Filtering::SaveFile(const char* ResultPath, unsigned char* Output)
{
	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//定义波段排列顺序
	int Bandmap[1] = { 1 };

	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* saveDataset = Driver->Create(ResultPath, this->Xsize, this->Ysize, this->Bandnum, this->DataType, papszOptions);

	//若创建失败
	if (!saveDataset)
	{
		cout << "create copy failed!" << endl;
	}
	else
		cout << "生成成功，请查看文件！" << endl;
	//使用RasterIO函数将数据写到创建的数据中
	saveDataset->RasterIO(GF_Write, 0, 0, this->Xsize, this->Ysize, Output, this->Xsize, this->Ysize, this->DataType, this->Bandnum, Bandmap, 0, 0, 0);

	GDALClose(saveDataset);
}