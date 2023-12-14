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
	std::cout << "�����������" << std::endl;
}

float min(float a, float b)
{
	return a <= b ? a : b;
}

void Tttt(const char* path1, const char* path2)
{
	//�������ú�����·������
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	//��ȡӰ�����ݼ�
	GDALDataset* A = (GDALDataset*)GDALOpen(path1, GA_ReadOnly); //����ף����� RGB ת���� HSI
	GDALDataset* B = (GDALDataset*)GDALOpen(path2, GA_ReadOnly); //ȫɫ��ֻ��һ�����Σ������滻 I
	//��ȡӰ����ߡ�����������������
	int A_Width = A->GetRasterXSize();
	int A_Height = A->GetRasterYSize();
	int A_BandCount = A->GetRasterCount();
	GDALDataType A_Type = A->GetRasterBand(1)->GetRasterDataType();

	int B_Width = B->GetRasterXSize();
	int B_Height = B->GetRasterYSize();
	int B_BandCount = B->GetRasterCount();
	GDALDataType B_Type = B->GetRasterBand(1)->GetRasterDataType();

	//�����ڴ�
	unsigned short* A_DataR = new unsigned short[A_Width * A_Height * A_Type];
	unsigned short* A_DataG = new unsigned short[A_Width * A_Height * A_Type];
	unsigned short* A_DataB = new unsigned short[A_Width * A_Height * A_Type];


	GDALRasterBand* aBandMulR;
	GDALRasterBand* aBandMulG;
	GDALRasterBand* aBandMulB;
	GDALRasterBand* bBandMul;
	double* B_Data = new double[B_Width * B_Height * B_Type];


	//��ȡ��Ӧ��������
	aBandMulB = A->GetRasterBand(1);
	aBandMulG = A->GetRasterBand(2);
	aBandMulR = A->GetRasterBand(3);
	aBandMulR->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataR, A_Width, A_Height, A_Type, 0, 0);
	aBandMulG->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataG, A_Width, A_Height, A_Type, 0, 0);
	aBandMulB->RasterIO(GF_Read, 0, 0, A_Width, A_Height, A_DataB, A_Width, A_Height, A_Type, 0, 0); 

	bBandMul = B->GetRasterBand(1);
	bBandMul->RasterIO(GF_Read, 0, 0, B_Width, B_Height, B_Data, B_Width, B_Height, B_Type, 0, 0);

	//�������ڱ��� H��S��I���ڴ�����
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

	//����Ӱ��ʹ�ֱ��ʴﵽ4096
	//ʹ��˫�����ڲ巨
	for (int y = 0; y < B_Height; y++)
	{
		for (int x = 0; x < B_Width; x++)
		{
			NewI_Data[y * B_Width + x] = B_Data[y * B_Width + x];

			float tmpX = x * (A_Width / B_Width);
			float tmpY = y * (A_Height / B_Height);
			//floor-����ȡ����ceil-����ȡ��
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


	//��任 HSI����>RGB
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




	//����Ӱ��ĵ�ַ
	const char* ResultPath = ".\\data\\result4.tif";

	//��ȡ����
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");

	//���岨������˳��
	int Bandmap[3] = { 1,2,3 };
	//��������Ӱ������ݼ�
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* outDataset = Driver->Create(ResultPath, B_Width, B_Height, 3, A_Type, papszOptions);

	//������ʧ��
	if (!outDataset)
	{
		std::cout << "Create File Failed!" << std::endl;
	}
	else
	{
		//ʹ��RasterIO����������д��������������
		//����B-G-R��˳��
		outDataset->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewB_Data, B_Width, B_Height, B_Type, 0, 0);
		outDataset->GetRasterBand(2)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewG_Data, B_Width, B_Height, B_Type, 0, 0);
		outDataset->GetRasterBand(3)->RasterIO(GF_Write, 0, 0, B_Width, B_Height, NewR_Data, B_Width, B_Height, B_Type, 0, 0);

	}

	//ɾ��������ڴ�
	delete[]A_DataR;

	

}

