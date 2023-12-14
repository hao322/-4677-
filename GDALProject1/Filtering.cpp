#include "pch.h"
#include "Filtering.h"
#include <gdal_priv.h>
#include <iostream>
#include <algorithm>

using namespace std;

//�вι��캯������ֵ
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

//��ֵ�˲�
unsigned char* Filtering::Mean()
{
	//�����ڴ� һ������
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

	int temp[9] = { 0 };//3*3��
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

//��ֵ�˲�
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
	int temp[9] = { 0 };//3*3��
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

//����Ӱ��
void Filtering::SaveFile(const char* ResultPath, unsigned char* Output)
{
	//��ȡ����
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//���岨������˳��
	int Bandmap[1] = { 1 };

	//��������Ӱ������ݼ�
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* saveDataset = Driver->Create(ResultPath, this->Xsize, this->Ysize, this->Bandnum, this->DataType, papszOptions);

	//������ʧ��
	if (!saveDataset)
	{
		cout << "create copy failed!" << endl;
	}
	else
		cout << "���ɳɹ�����鿴�ļ���" << endl;
	//ʹ��RasterIO����������д��������������
	saveDataset->RasterIO(GF_Write, 0, 0, this->Xsize, this->Ysize, Output, this->Xsize, this->Ysize, this->DataType, this->Bandnum, Bandmap, 0, 0, 0);

	GDALClose(saveDataset);
}