#include "pch.h"
#include "HIS_ImgFusion.h"
#include <iostream>
const double pi = 4 * atan(1);


HIS_ImgFusion::HIS_ImgFusion(const char* DataPath)
{
	this->Dataset = (GDALDataset*)GDALOpen(DataPath, GA_ReadOnly);
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
	this->DataType = Dataset->GetRasterBand(1)->GetRasterDataType();
}
int HIS_ImgFusion::getXsize()
{
	return this->Xsize;
}
int HIS_ImgFusion::getYsize()
{
	return this->Ysize;
}
int HIS_ImgFusion::getBandnum()
{
	return this->Bandnum;
}
HIS_ImgFusion::~HIS_ImgFusion()
{
	std::cout << "程序清理完成" << std::endl;
}

float min(float a, float b)
{
	return a <= b ? a : b;
}

void Tttt(const char* path1, const char* path2)
{
	//驱动设置和中文路径设置
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	//读取影像数据集
	GDALDataset* A = (GDALDataset*)GDALOpen(path1, GA_ReadOnly); //多光谱，用于 RGB 转换到 HSI
	GDALDataset* B = (GDALDataset*)GDALOpen(path2, GA_ReadOnly); //全色，只有一个波段，用于替换 I
	//获取影像宽、高、波段数及数据类型
	int A_Width = A->GetRasterXSize();
	int A_Height = A->GetRasterYSize();
	int A_BandCount = A->GetRasterCount();
	GDALDataType A_Type = A->GetRasterBand(1)->GetRasterDataType();

	int B_Width = B->GetRasterXSize();
	int B_Height = B->GetRasterYSize();
	int B_BandCount = B->GetRasterCount();
	GDALDataType B_Type = B->GetRasterBand(1)->GetRasterDataType();

	//申请内存
	unsigned short* A_DataR = new unsigned short[A_Width * A_Height * A_Type];
	unsigned short* A_DataG = new unsigned short[A_Width * A_Height * A_Type];
	unsigned short* A_DataB = new unsigned short[A_Width * A_Height * A_Type];


	GDALRasterBand* aBandMulR;
	GDALRasterBand* aBandMulG;
	GDALRasterBand* aBandMulB;
	GDALRasterBand* bBandMul;
	double* B_Data = new double[B_Width * B_Height * B_Type];


	//读取对应波段数据
	aBandMulB = A->GetRasterBand(1);
	aBandMulG = A->GetRasterBand(2);
	aBandMulR = A->GetRasterBand(3);
	aBandMulR->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataR, A_Width, A_Height, A_Type, 0, 0);
	aBandMulG->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataG, A_Width, A_Height, A_Type, 0, 0);
	aBandMulB->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataB, A_Width, A_Height, A_Type, 0, 0); 

	bBandMul = B->GetRasterBand(1);
	bBandMul->RasterIO(GF_Read, 0, 0, B_Width, B_Height, B_Data, B_Width, B_Height, B_Type, 0, 0);

	//申请用于保存 H、S、I的内存数组
	/*unsigned float* theta = new unsigned float[A_Width * A_Height * A_Type];
	unsigned float* H_Data = new unsigned float[A_Width * A_Height * A_Type];
	unsigned float* S_Data = new unsigned float[A_Width * A_Height * A_Type];*/
	float* theta = new float[A_Width * A_Height * A_Type];
	float* H_Data = new float[A_Width * A_Height * A_Type];
	float* S_Data = new float[A_Width * A_Height * A_Type];
	float* I_Data = new float[A_Width * A_Height * A_Type];

	//unsigned short* I_Data = new unsigned short[A_Width * A_Height * A_Type];
	/*unsigned short* NewI_Data = new unsigned short[B_Width * B_Height * B_Type];
	unsigned short* NewH_Data = new unsigned short[B_Width * B_Height * B_Type];
	unsigned short* NewS_Data = new unsigned short[B_Width * B_Height * B_Type];*/
	float* NewH_Data = new float[A_Width * A_Height * A_Type];
	float* NewS_Data = new float[A_Width * A_Height * A_Type];
	float* NewI_Data = new float[A_Width * A_Height * A_Type];

	unsigned short* NewR_Data = new unsigned short[B_Width * B_Height];
	unsigned short* NewG_Data = new unsigned short[B_Width * B_Height];
	unsigned short* NewB_Data = new unsigned short[B_Width * B_Height];




	int k = 0;
	for (int i = 0; i < A_Height; i++) // Height = Ysize
	{
		for (int j = 0; j < A_Width; j++) // Width = Xsize
		{
			k = i * A_Width + j;

			S_Data[k] = 1 - 3 * min(min(A_DataR[k], A_DataB[k]), A_DataG[k]) / (A_DataR[k] + A_DataB[k] + A_DataG[k]);
			I_Data[k] = (A_DataR[k] + A_DataB[k] + A_DataG[k]) / 3;

			theta[k] = acos(0.5 * (2 * A_DataR[k] - A_DataG[k] - A_DataB[k]) / 
				sqrt(pow(A_DataR[k] - A_DataG[k], 2) + (A_DataR[k] - A_DataB[k]) * (A_DataG[k] - A_DataB[k])));
			theta[k] = theta[k] * 180 / pi;

			if (A_DataB[k] <= A_DataG[k])
			{
				H_Data[k] = theta[k];
			}
			else
			{
				H_Data[k] = 360 - theta[k];
			}
		}
	}

	//拉伸影像，使分辨率达到4096
	//使用双线性内插法
	for (int y = 0; y < B_Height; y++)
	{
		for (int x = 0; x < B_Width; x++)
		{
			NewI_Data[y * B_Width + x] = B_Data[y * B_Width + x];

			float tmpX = x * (A_Width / B_Width);
			float tmpY = y * (A_Height / B_Height);
			//floor-向下取整，ceil-向上取整
			int flr_X = floor(tmpX), flr_Y = floor(tmpY);
			int cel_X = ceil(tmpX), cel_Y = ceil(tmpY);
			float d_t = tmpY - floor(tmpY);
			float d_b = ceil(tmpY) - tmpY;
			float d_l = tmpX - floor(tmpX);
			float d_r = ceil(tmpX) - tmpX;

			NewH_Data[y * B_Width + x] = H_Data[flr_Y * B_Width + flr_X] * d_b * d_r + H_Data[flr_Y * B_Width + cel_X] * d_b * d_l +
				H_Data[cel_Y * B_Width + flr_X] * d_t * d_r + H_Data[cel_Y * B_Width + cel_X] * d_t * d_l;
			NewS_Data[y * B_Width + x] = S_Data[flr_Y * B_Width + flr_X] * d_b * d_r + S_Data[flr_Y * B_Width + cel_X] * d_b * d_l +
				S_Data[cel_Y * B_Width + flr_X] * d_t * d_r + S_Data[cel_Y * B_Width + cel_X] * d_t * d_l;
		}
	}


	//逆变换 HSI——>RGB
	for (int i = 0; i < B_Height; i++)
	{
		for (int j = 0; j < B_Width; j++)
		{
			k = i * B_Width + j;
			float h = NewH_Data[k];
			float H = NewH_Data[k] * pi /180 ;
			float S = NewS_Data[k];
			float I = NewI_Data[k];

			if (h >= 0 && h < 120)
			{
				NewB_Data[k] = I * (1 - S);
				NewR_Data[k] = I * (1 + S * cos(H) / cos(pi / 3 - H));
				NewG_Data[k] = 3 * I - (NewB_Data[k] + NewR_Data[k]);
				
			}
			else if ( h >= 120 && H < 240 )
			{
				H = H - pi * 2 / 3;
				NewR_Data[k] = I * (1 - S);
				NewG_Data[k] = I * (1 + S * cos(H) / cos(pi / 3 - H));
				NewB_Data[k] = 3 * I - (NewR_Data[k] + NewG_Data[k]);
			}
			else
			{
				H = H - pi * 4 / 3;
				NewG_Data[k] = I * (1 - S);
				NewB_Data[k] = I * (1 + S * cos(H) / cos(pi / 3 - H));
			}
			NewB_Data[k] = NewB_Data[k] * 255;
			NewR_Data[k] = NewR_Data[k] * 255;
			NewG_Data[k] = NewG_Data[k] * 255;
		}
	}




	//保存影像的地址
	const char* ResultPath = ".\\data\\result4.tif";

	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	//定义波段排列顺序
	int Bandmap[3] = { 1,2,3 };
	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* outDataset = Driver->Create(ResultPath, B_Width, B_Height, 3, A_Type, papszOptions);

	//若创建失败
	if (!outDataset)
	{
		std::cout << "Create File Failed!" << std::endl;
	}
	else
	{
		//使用RasterIO函数将数据写到创建的数据中
		//按照B-G-R的顺序
		outDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewB_Data, B_Width, B_Height, B_Type, 0, 0);
		outDataset->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewG_Data, B_Width, B_Height, B_Type, 0, 0);
		outDataset->GetRasterBand(3)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewR_Data, B_Width, B_Height, B_Type, 0, 0);

	}

	//删除申请的内存
	delete[]A_DataR;

	

}

