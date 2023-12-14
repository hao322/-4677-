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
	//�����ڴ棬��������
	unsigned short* Nir_Img = new unsigned short[Xsize * Ysize * DataType];
	unsigned short* Red_Img = new unsigned short[Xsize * Ysize * DataType]; 
	 
	//�����ڴ棬���ڱ���NDVI
	float* NDVI = new float[Xsize * Ysize * sizeof(float)];

	//������
	GDALRasterBand* band_NIR = Dataset->GetRasterBand(4);
	band_NIR->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Nir_Img, Xsize, Ysize, DataType, 0, 0);
	//�첨��
	GDALRasterBand* band_Red = Dataset->GetRasterBand(3);
	band_Red->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Red_Img, Xsize, Ysize, DataType, 0, 0);

	//����Ӱ�񣬼���NDVI
	for (int i = 0; i < Ysize; i++)   
	{
		for (int j = 0; j < Xsize; j++)
		{
			//��ĸ����һ����С���������ĸΪ0
			NDVI[Xsize * i + j] = float((Nir_Img[Xsize * i + j] - Red_Img[Xsize * i + j]) / (1.0e-7 + Nir_Img[Xsize * i + j] + Red_Img[Xsize * i + j]));

			//��ֵ�ָ����ֵ��
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
	//��ȡ����
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//���岨������˳��
	int Bandmap[1] = { 1 };

	//��������Ӱ������ݼ�
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* saveDataset = Driver->Create(ResultPath, this->Xsize, this->Ysize, 1, GDT_Float32, papszOptions);

	//������ʧ��
	if (!saveDataset)
	{
		cout << "Create file failed!" << endl;
	}
	else
	{
		//ʹ��RasterIO����������д�������������У�ע���ʱֻ��һ������
		saveDataset->RasterIO(GF_Write, 0, 0, this->Xsize, this->Ysize, NDVI, this->Xsize, this->Ysize, GDT_Float32, 1, Bandmap, 0, 0, 0);

		cout << "Success����鿴�ļ���" << endl;
	}

	GDALClose(Dataset);
	GDALClose(saveDataset);
}