#include "pch.h"
#include "Data.h"
#include <iostream>

using namespace std;

Data::Data(const char* DataPath) 
{
	this->Dataset = (GDALDataset*)GDALOpen(DataPath, GA_ReadOnly);
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
	this->DataType = Dataset->GetRasterBand(1)->GetRasterDataType();
}

int Data::getXsize()
{
	return this->Xsize;
}
int Data::getBandnum()
{
	return this->Bandnum;
}

float* Data::getNDVI()
{
	//申请内存，两个波段
	unsigned short* Nir_Img = new unsigned short[Xsize * Ysize * DataType];
	unsigned short* Red_Img = new unsigned short[Xsize * Ysize * DataType]; 
	 
	//申请内存，用于保存NDVI
	float* NDVI = new float[Xsize * Ysize * sizeof(float)];

	//近红外
	GDALRasterBand* band_NIR = Dataset->GetRasterBand(4);
	band_NIR->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Nir_Img, Xsize, Ysize, DataType, 0, 0);
	//红波段
	GDALRasterBand* band_Red = Dataset->GetRasterBand(3);
	band_Red->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Red_Img, Xsize, Ysize, DataType, 0, 0);

	//遍历影像，计算NDVI
	for (int i = 0; i < Ysize; i++)   
	{
		for (int j = 0; j < Xsize; j++)
		{
			//分母加上一个极小量，避免分母为0
			NDVI[Xsize * i + j] = float((Nir_Img[Xsize * i + j] - Red_Img[Xsize * i + j]) / (1.0e-7 + Nir_Img[Xsize * i + j] + Red_Img[Xsize * i + j]));

			//阈值分割——二值化
			if (NDVI[Xsize * i + j] > 0.3)
			{
				NDVI[Xsize * i + j] = 1;
			}
			else
			{
				NDVI[Xsize * i + j] = -1;
			}

		}
	}
	return NDVI;
}

void Data::saveFile(const char* ResultPath, float* NDVI)
{
	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//定义波段排列顺序
	int Bandmap[1] = { 1 };

	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* saveDataset = Driver->Create(ResultPath, this->Xsize, this->Ysize, 1, GDT_Float32, papszOptions);

	//若创建失败
	if (!saveDataset)
	{
		cout << "Create file failed!" << endl;
	}
	else
	{
		//使用RasterIO函数将数据写到创建的数据中，注意此时只有一个波段
		saveDataset->RasterIO(GF_Write, 0, 0, this->Xsize, this->Ysize, NDVI, this->Xsize, this->Ysize, GDT_Float32, 1, Bandmap, 0, 0, 0);

		cout << "Success，请查看文件！" << endl;
	}

	GDALClose(Dataset);
	GDALClose(saveDataset);
}