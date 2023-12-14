#include "pch.h"
#include "Data.h"
#include <gdal_priv.h>

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
int Data::getYsize()
{
	return this->Ysize;
}
int Data::getBandnum()
{
	return this->Bandnum;
}

pair<int, int> Data::match_ncc(Data &tem)
{
	//匹配位置、ncc
	int r, c;
	float ncc = -1;

	//申请内存，存储两张图像的波段
	unsigned char* Img = new unsigned char[Xsize * Ysize * DataType];
	unsigned char* Tem = new unsigned char[tem.Xsize * tem.Ysize * tem.DataType];
	//用于存储模板影像在待匹配影像滑动时的小影像数组
	unsigned char* buffer = new unsigned char[tem.Xsize * tem.Ysize * tem.DataType];

	//图像IMG
	GDALRasterBand* band_Img = Dataset->GetRasterBand(1);
	band_Img->RasterIO(GF_Read, 0, 0, Xsize, Ysize, Img, Xsize, Ysize, DataType, 0, 0);

	//模板影像
	GDALRasterBand* band_Tem = tem.Dataset->GetRasterBand(1);
	band_Tem->RasterIO(GF_Read, 0, 0, tem.Xsize, tem.Ysize, Tem, tem.Xsize, tem.Ysize, tem.DataType, 0, 0);


	//Calculate
	for (int i = 0; i < Ysize - tem.Ysize; i++)
	{
		for (int j = 0; j < Xsize - tem.Xsize; j++)
		{
			for (int k = 0; k < tem.Xsize * tem.Ysize; k++)
			{
				*(buffer + k) = Img[Xsize * (i + k / tem.Xsize) + j + k % tem.Xsize];
			}
			float temp = cal_NCC(buffer, Tem, tem.Xsize * tem.Ysize);
			if (temp > ncc)
			{
				ncc = temp;
				r = i;
				c = j;
			}

		}

	}
	return make_pair(r, c);
}

float Data::cal_Mean(unsigned char* array, int length)
{
	float sum = 0;
	for (int i = 0; i < length; i++)
	{
		sum += *(array + i);
	}
	return sum / length;
}

float Data::cal_NCC(unsigned char* tem, unsigned char* img, int size)
{
	float tem_Mean = cal_Mean(tem, sizeof(tem));
	float img_Mean = cal_Mean(img, sizeof(img));

	float sum_n1 = 0;
	float sum_n21 = 0;
	float sum_n22 = 0;

	for (int i = 0; i < size; i++)
	{
		sum_n1 += (*(tem + i) - tem_Mean) * (*(img + i) - img_Mean);
		sum_n21 += pow((*(tem + i) - tem_Mean), 2);
		sum_n22 += pow((*(img + i) - img_Mean), 2);
	}

	return sum_n1 / sqrt(sum_n21 * sum_n22);
}


