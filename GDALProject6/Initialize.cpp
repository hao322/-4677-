#include "pch.h"
#include "Initialize.h"

 
Initialize::Initialize(const char* datapath) 
{
	this->Dataset = (GDALDataset*)GDALOpen(datapath, GA_ReadOnly);
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
	this->Datatype = Dataset->GetRasterBand(1)->GetRasterDataType();
}
int Initialize::getXsize()
{
	return this->Xsize;
}
int Initialize::getYsize()
{
	return this->Ysize;
}
int Initialize::getBandnum()
{
	return this->Bandnum;
}
GDALDataset* Initialize::getDataset()
{
	return this->Dataset;
}
GDALDataType Initialize::getDatatype()
{
	return this->Datatype;
}
Initialize::~Initialize() {};

unsigned int* getChange(Initialize& after, Initialize& before)
{
	//after
	GDALDataType a_type = after.getDatatype(); //GDT_UInt16 对应 unsigned short
	GDALDataset* a_Dataset = after.getDataset();
	int a_Xsize = after.getXsize();
	int a_Ysize = after.getYsize();
	int a_Bandnum = after.getBandnum();
	//before
	GDALDataType b_type = before.getDatatype();
	GDALDataset* b_Dataset = before.getDataset();
	int b_Xsize = before.getXsize();
	int b_Ysize = before.getYsize();
	int b_Bandnum = before.getBandnum();
	
	unsigned int* mChange = new unsigned int[a_Xsize * a_Ysize];

	if (a_Bandnum != b_Bandnum || a_Xsize != b_Xsize || a_Ysize != b_Ysize)
	{
		std::cout << "两图像不匹配，无法进行变化检测" << std::endl;
	}
	else
	{
		GDALRasterBand* befBand1 = b_Dataset->GetRasterBand(1);
		GDALRasterBand* befBand2 = b_Dataset->GetRasterBand(2);
		GDALRasterBand* befBand3 = b_Dataset->GetRasterBand(3);
		GDALRasterBand* befBand4 = b_Dataset->GetRasterBand(4);

		GDALRasterBand* aftBand1 = a_Dataset->GetRasterBand(1);
		GDALRasterBand* aftBand2 = a_Dataset->GetRasterBand(2);
		GDALRasterBand* aftBand3 = a_Dataset->GetRasterBand(3);
		GDALRasterBand* aftBand4 = a_Dataset->GetRasterBand(4);

		unsigned short* aBand1 = new unsigned short[a_Xsize * a_Ysize];
		unsigned short* aBand2 = new unsigned short[a_Xsize * a_Ysize];
		unsigned short* aBand3 = new unsigned short[a_Xsize * a_Ysize];
		unsigned short* aBand4 = new unsigned short[a_Xsize * a_Ysize];

		unsigned short* bBand1 = new unsigned short[b_Xsize * b_Ysize];
		unsigned short* bBand2 = new unsigned short[b_Xsize * b_Ysize];
		unsigned short* bBand3 = new unsigned short[b_Xsize * b_Ysize];
		unsigned short* bBand4 = new unsigned short[b_Xsize * b_Ysize];

		befBand1->RasterIO(GF_Read, 0, 0, b_Xsize, b_Ysize, bBand1, b_Xsize, b_Ysize, b_type, 0, 0);
		befBand2->RasterIO(GF_Read, 0, 0, b_Xsize, b_Ysize, bBand2, b_Xsize, b_Ysize, b_type, 0, 0);
		befBand3->RasterIO(GF_Read, 0, 0, b_Xsize, b_Ysize, bBand3, b_Xsize, b_Ysize, b_type, 0, 0);
		befBand4->RasterIO(GF_Read, 0, 0, b_Xsize, b_Ysize, bBand4, b_Xsize, b_Ysize, b_type, 0, 0);

		aftBand1->RasterIO(GF_Read, 0, 0, a_Xsize, a_Ysize, aBand1, a_Xsize, a_Ysize, a_type, 0, 0);
		aftBand2->RasterIO(GF_Read, 0, 0, a_Xsize, a_Ysize, aBand2, a_Xsize, a_Ysize, a_type, 0, 0);
		aftBand3->RasterIO(GF_Read, 0, 0, a_Xsize, a_Ysize, aBand3, a_Xsize, a_Ysize, a_type, 0, 0);
		aftBand4->RasterIO(GF_Read, 0, 0, a_Xsize, a_Ysize, aBand4, a_Xsize, a_Ysize, a_type, 0, 0);

		int mark = 0;
		for (int i = 0; i < a_Ysize; i++)
		{
			for (int j = 0; j < a_Xsize; j++)
			{
				mark = i * a_Xsize + j;

				double temp1 = (double)aBand1[mark] - (double)bBand1[mark];
				double temp2 = (double)aBand2[mark] - (double)bBand2[mark];
				double temp3 = (double)aBand3[mark] - (double)bBand3[mark];
				double temp4 = (double)aBand4[mark] - (double)bBand4[mark];

				double value = sqrt(temp1 * temp1 + temp2 * temp2 + temp3 * temp3 + temp4 * temp4);

				/*if (value > 255)
				{
					mChange[mark] = 255;
				}
				else
				{
					mChange[mark] = int(value + 0.5);
				}*/
			}
		}
	}
	return mChange;
}
//确定阈值
int threshold(unsigned int* change, int Xsize, int Ysize)
{
	int Bin[256];
	double Percent[256];

	for (int i = 0; i < 256; i++)
	{
		int cnt = 0;
		for (int j = 0; j < Xsize * Ysize; j++)
		{
			if (change[j] == i)
			{
				cnt++;
			}

			Bin[i] = cnt;
			Percent[i] = cnt / (Xsize * Ysize);
		}
	}

	int T = 0;
	double G = 0;

	for (int t = 0; t < 256; t++)
	{
		double N0 = 0.0, N1 = 0.0;
		double sum0 = 0, sum1 = 0;

		for (int i = 0; i < Xsize * Ysize; i++)
		{
			if (change[i] < t)
			{
				N0++;
				sum0 += change[i];
			}
			else
			{
				N1++;
				sum1 += change[i];
			}
		}

		double u0 = sum0 / N0;
		double u1 = sum1 / N1;

		double w0 = N0 / (Xsize * Ysize);
		double w1 = N1 / (Xsize * Ysize);

		double g = w0 * w1 * (u0 - u1) * (u0 - u1);

		if (g > G)
		{
			G = g;
			T = t;
		}
	}

	return T;

}

//二值化
unsigned short* segmentation(unsigned int* change, int threshold, int Xsize, int Ysize)
{
	unsigned short* result = new unsigned short[Xsize * Ysize];
	for (int i = 0; i < Ysize; i++)
	{
		for (int j = 0; j < Xsize; j++)
		{
			int mark = i * Xsize + j;
			if (change[mark] <= threshold)
			{
				result[mark] = 0;
			}
			else
			{
				result[mark] = 255;
			}

		}
	}
	return result;
}

void saveFile(unsigned short* result, const char* resultPath, Initialize aim)
{
	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//定义波段排列顺序
	int Bandmap[1] = { 1 };

	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* saveDataset = Driver->Create(resultPath, aim.getXsize(), aim.getYsize(), 1, aim.getDatatype(), papszOptions);

	//若创建失败
	if (!saveDataset)
	{
		std::cout << "Create file failed!" << std::endl;
	}
	else
	{
		//使用RasterIO函数将数据写到创建的数据中
		saveDataset->RasterIO(GF_Write, 0, 0, aim.getXsize(), aim.getYsize(), result, aim.getXsize(), aim.getYsize(), aim.getDatatype(), 1, Bandmap, 0, 0, 0);

		std::cout << "Success，请查看文件！" << std::endl;
	}

	GDALClose(saveDataset);
}


//精度评定
void EvaAccuracy(Initialize& ref, Initialize& res)
{
	//reference
	GDALDataType ref_type = ref.getDatatype(); //GDT_UInt16 对应 unsigned short
	GDALDataset* ref_Dataset = ref.getDataset();
	int ref_Xsize = ref.getXsize();
	int ref_Ysize = ref.getYsize();
	int ref_Bandnum = ref.getBandnum();
	//result
	GDALDataType res_type = res.getDatatype();
	GDALDataset* res_Dataset = res.getDataset();
	int res_Xsize = res.getXsize();
	int res_Ysize = res.getYsize();
	int res_Bandnum = res.getBandnum();

	/*if (ref_type == GDT_UInt16)
	{
		unsigned short* _resBand = new unsigned short[ref_Xsize * ref_Ysize];
	}
	else if (ref_type == GDT_Byte)
	{
		unsigned char* _resBand = new unsigned char[ref_Xsize * ref_Ysize];

	}*/

	GDALRasterBand* refBand = ref_Dataset->GetRasterBand(1);
	unsigned short* _refBand = new unsigned short[ref_Xsize * ref_Ysize * ref_type];
	refBand->RasterIO(GF_Read, 0, 0, ref_Xsize, ref_Ysize, _refBand, ref_Xsize, ref_Ysize, ref_type, 0, 0);

	GDALRasterBand* resBand = res_Dataset->GetRasterBand(1);
	unsigned short* _resBand = new unsigned short[res_Xsize * res_Ysize * res_type];
	resBand->RasterIO(GF_Read, 0, 0, res_Xsize, res_Ysize, _resBand, res_Xsize, res_Ysize, res_type, 0, 0);


	int N11 = 0; int N12 = 0; int N21 = 0; int N22 = 0;

	int mark = 0;
	for (int i = 0; i < ref_Xsize * ref_Ysize; i++)
	{
		if (_refBand[i] == 255 && _resBand[i] == 255)
		{
			N11++;
		}
		else if (_refBand[i] == 255 && _resBand[i] == 0)
		{
			N21++;
		}
		else if (_refBand[i] == 0 && _resBand[i] == 255)
		{
			N12++;
		}
		else
		{
			N22++;
		}
	}

	//混淆矩阵
	int A_refChange = N11 + N21;
	int B_refNoChange = N12 + N22;
	int C_resChange = N11 + N12;
	int D_resNoChange = N21 + N22;
	int T_Count = A_refChange + B_refNoChange + C_resChange + D_resNoChange;

	//计算漏检率、虚警率、总分类精度、Kappa系数
	double p1 = N21 / A_refChange;
	double p2 = N12 / C_resChange;
	double p = (N11 + N22) / T_Count;
	double kappa = (T_Count * (N11 + N22) - (A_refChange * C_resChange + B_refNoChange * D_resNoChange)) /
		(T_Count * T_Count - (A_refChange * C_resChange + B_refNoChange * D_resNoChange));


}