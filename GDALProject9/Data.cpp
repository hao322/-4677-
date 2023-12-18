#include "pch.h"
#include "Data.h"

// ͼ��Ŀ�Ⱥ͸߶�
const int Width = 500;
const int Height = 1000;

const double PI = 4 * atan(1);


Data::Data(const char* FilePath)
{
	this->Dataset = (GDALDataset*)GDALOpen(FilePath, GA_ReadOnly);
	this->Datatype = Dataset->GetRasterBand(1)->GetRasterDataType();
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
}
int Data::getXsize()
{
	return Xsize;
}
int Data::getYsize()
{
	return Ysize;
}
int Data::getBandnum()
{
	return Bandnum;
}
GDALDataset* Data::getDataset()
{
	return Dataset;
}
GDALDataType Data::getDatatype()
{
	return Datatype;
}
Data::~Data() {};

// ����������������
double caldistance(Pixel a, Pixel b)
{
	return sqrt(pow(a.r - b.r, 2) + pow(a.g - b.g, 2) + pow(a.b - b.b, 2));
}
// ��ȡͼ������
vector<Pixel> loadImage(Data& Img)
{
	using namespace std;
	vector<Pixel> pl;
	GDALDataType ImgType = Img.getDatatype();
	// ��������
	GDALDataset* imgDataSet = Img.getDataset();
	unsigned char* PixelRed = new unsigned char[Width * Height];
	unsigned char* PixelGreen = new unsigned char[Width * Height];
	unsigned char* PixelBlue = new unsigned char[Width * Height];

	GDALRasterBand* Bandr = imgDataSet->GetRasterBand(1);
	GDALRasterBand* Bandg = imgDataSet->GetRasterBand(2);
	GDALRasterBand* Bandb = imgDataSet->GetRasterBand(3);

	Bandr->RasterIO(GF_Read, 0, 0, Width, Height, PixelRed, Width, Height, ImgType, 0, 0);
	Bandg->RasterIO(GF_Read, 0, 0, Width, Height, PixelGreen, Width, Height, ImgType, 0, 0);
	Bandb->RasterIO(GF_Read, 0, 0, Width, Height, PixelBlue, Width, Height, ImgType, 0, 0);

	int Mark = 0;
	for (int i = 0; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			Mark = i * Width + j;
			Pixel pix;
			pix.r = PixelRed[Mark];
			pix.g = PixelGreen[Mark];
			pix.b = PixelBlue[Mark];
			pl.push_back(pix);
		}
	}
	delete[] PixelRed;
	delete[] PixelGreen;
	delete[] PixelBlue;
	cout << "load Image ... over" << endl;
	return pl;
}

// ������ֵ��Э����������
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeans(vector<Pixel> train, vector<Pixel> Label)
{
	vector<Pixel> meanPixels(5);
	int mark = 0; int k0 = 0; int k1 = 0; int k2 = 0; int k3 = 0; int k4 = 0;

	double sumR0 = 0;	double sumG0 = 0;	double sumB0 = 0;
	double sumR50 = 0;	double sumG50 = 0;	double sumB50 = 0;
	double sumR100 = 0;	double sumG100 = 0;	double sumB100 = 0;
	double sumR150 = 0;	double sumG150 = 0;	double sumB150 = 0;
	double sumR200 = 0;	double sumG200 = 0;	double sumB200 = 0;

	for (int i = 0; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			mark = i * Width + j;
			if(Label[mark].r == 0)
			{
				sumR0 += train[mark].r;
				sumG0 += train[mark].g;
				sumB0 += train[mark].b;
				k0++;
			}
			else if(Label[mark].r == 50)
			{
				sumR50 += train[mark].r;
				sumG50 += train[mark].g;
				sumB50 += train[mark].b;
				k1++;
			}
			else if (Label[mark].r == 100)
			{
				sumR100 += train[mark].r;
				sumG100 += train[mark].g;
				sumB100 += train[mark].b;
				k2++;
			}
			else if (Label[mark].r == 150)
			{
				sumR150 += train[mark].r;
				sumG150 += train[mark].g;
				sumB150 += train[mark].b;
				k3++;
			}
			else
			{
				sumR200 += train[mark].r;
				sumG200 += train[mark].g;
				sumB200 += train[mark].b;
				k4++;
			}
		}
	}

	vector<double> probability(5);
	probability[0] = (double)k0 / train.size();
	probability[1] = (double)k1 / train.size();
	probability[2] = (double)k2 / train.size();
	probability[3] = (double)k3 / train.size();
	probability[4] = (double)k4 / train.size();

	meanPixels[0].r = sumR0 / k0; meanPixels[0].g = sumG0 / k0;	meanPixels[0].b = sumB0 / k0;
	meanPixels[1].r = sumR50 / k1; meanPixels[1].g = sumG50 / k1; meanPixels[1].b = sumB50 / k1;
	meanPixels[2].r = sumR100 / k2; meanPixels[2].g = sumG100 / k2;	meanPixels[2].b = sumB100 / k2;
	meanPixels[3].r = sumR150 / k3; meanPixels[3].g = sumG150 / k3;	meanPixels[3].b = sumB150 / k3;
	meanPixels[4].r = sumR200 / k4;  meanPixels[4].g = sumG200 / k4; meanPixels[4].b = sumB200 / k4;

	vector<Vector3d> means;
	for (int i = 0; i < meanPixels.size(); i++)
	{
		Vector3d temp = Pix2Vec(meanPixels[i]);
		means.push_back(temp);
	}

	MatrixXd cov0(k0, 3);	MatrixXd cov1(k1, 3);	MatrixXd cov2(k2, 3);	MatrixXd cov3(k3, 3);	MatrixXd cov4(k4, 3);
	int a0 = 0; int a1 = 0; int a2 = 0; int a3 = 0; int a4 = 0;

	for (int i = 0; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			mark = i * Width + j;
			if (Label[mark].r == 0)
			{
				cov0(a0, 0) = train[mark].r - meanPixels[0].r;
				cov0(a0, 1) = train[mark].g - meanPixels[0].g;
				cov0(a0, 2) = train[mark].b - meanPixels[0].b;
				a0++;
			}
			else if (Label[mark].r == 50)
			{
				cov1(a1, 0) = train[mark].r - meanPixels[1].r;
				cov1(a1, 1) = train[mark].g - meanPixels[1].g;
				cov1(a1, 2) = train[mark].b - meanPixels[1].b;
				a1++;
			}
			else if (Label[mark].r == 100)
			{
				cov2(a2, 0) = train[mark].r - meanPixels[2].r;
				cout << cov2(a2, 0) << endl;
				cov2(a2, 1) = train[mark].g - meanPixels[2].g;
				cov2(a2, 2) = train[mark].b - meanPixels[2].b;
				a2++;
			}
			else if (Label[mark].r == 150)
			{
				cov3(a3, 0) = train[mark].r - meanPixels[3].r;
				cov3(a3, 1) = train[mark].g - meanPixels[3].g;
				cov3(a3, 2) = train[mark].b - meanPixels[3].b;
				a3++;
			}
			else
			{
				cov4(a4, 0) = train[mark].r - meanPixels[4].r;
				cov4(a4, 1) = train[mark].g - meanPixels[4].g;
				cov4(a4, 2) = train[mark].b - meanPixels[4].b;
				a4++;
			}
		}
	}

	Matrix3d Cov0 = cov0.transpose() * cov0 * 0.5;
	Matrix3d Cov1 = cov1.transpose() * cov1 * 0.5;
	Matrix3d Cov2 = cov2.transpose() * cov2 * 0.5;
	Matrix3d Cov3 = cov3.transpose() * cov3 * 0.5;
	Matrix3d Cov4 = cov4.transpose() * cov4 * 0.5;
	
	vector<Matrix3d> Cov(5);
	Cov[0] = Cov0;
	Cov[1] = Cov1;
	Cov[2] = Cov2;
	Cov[3] = Cov3;
	Cov[4] = Cov4;

	// vector��0��4�ֱ��Ӧ0��50��100��150��200
	return make_tuple(means, Cov, probability);
}




// ��˹�����ܶȼ���
double gaussProDen(double means, double var, double x)
{
	double temp = -1 * (x - means) * (x - means) / (2 * var);
	double result = 1 / (sqrt(2 * PI * var)) * exp(temp);
	return result;
}

// ��Ԫ��˹�����ܶȼ���
double testgauss(VectorXd mean, MatrixXd cov, VectorXd test)
{
	int k = mean.size();

	double determinant = cov.determinant();
	double normalization = 1.0 / (pow(2 * M_PI, k / 2) * sqrt(determinant));

	MatrixXd inv_cov = cov.inverse();
	double exponent = -0.5 * (test - mean).transpose() * inv_cov * (test - mean);
	return normalization * exp(exponent);

}

// ���ص�������ת��
Vector3d Pix2Vec(Pixel a)
{
	VectorXd V(3);
	V << a.r, a.g, a.b;
	return V;
}

// �����Ȼ������
int testMax(vector<double> probbly, vector<Vector3d> means, vector<Matrix3d> cov, int trainLabels[], int numTrain, Pixel testPix)
{
	Vector3d tePix = Pix2Vec(testPix);

	int maxId = 0;
	double maxProb = 0;

	for (int i = 0; i < numTrain; i++)
	{
		double prob = testgauss(means[i], cov[i], tePix) * probbly[i];
		if (prob > maxProb)
		{
			maxProb = prob;
			maxId = trainLabels[i];
		}
	}
	return maxId;
}

// ��С�������
int minDistClassify(vector<Pixel> meanPix, int trainLabels[], int numTrain, Pixel testPix)
{
	double minDist = INT_MAX;
	unsigned int predLabel;

	for (int i = 0; i < numTrain; ++i)
	{
		// ����Ӧ�����������ֵ����֮��ľ���
		double dist = caldistance(testPix, meanPix[i]);

		if (dist < minDist)
		{
			minDist = dist;
			predLabel = trainLabels[i];
		}
	}
	// Ԥ�����
	return predLabel;
}

//��������
void EvaAccuracy(Data& ref, Data& res)
{
	//reference
	GDALDataType ref_type = ref.getDatatype(); //GDT_UInt16 ��Ӧ unsigned short
	GDALDataset* ref_Dataset = ref.getDataset();
	int ref_Xsize = ref.getXsize();
	int ref_Ysize = ref.getYsize();
	// int ref_Bandnum = ref.getBandnum();
	//result
	GDALDataType res_type = res.getDatatype();
	GDALDataset* res_Dataset = res.getDataset();
	int res_Xsize = res.getXsize();
	int res_Ysize = res.getYsize();
	// int res_Bandnum = res.getBandnum();

	GDALRasterBand* refBand = ref_Dataset->GetRasterBand(1);
	unsigned char* _refBand = new unsigned char[ref_Xsize * ref_Ysize];
	refBand->RasterIO(GF_Read, 0, 0, ref_Xsize, ref_Ysize, _refBand, ref_Xsize, ref_Ysize, ref_type, 0, 0);

	GDALRasterBand* resBand = res_Dataset->GetRasterBand(1);
	int* _resBand = new int[res_Xsize * res_Ysize];
	resBand->RasterIO(GF_Read, 0, 0, res_Xsize, res_Ysize, _resBand, res_Xsize, res_Ysize, res_type, 0, 0);


	int N11 = 0; int N12 = 0; int N13 = 0; int N14 = 0; int N15 = 0;
	int N21 = 0; int N22 = 0; int N23 = 0; int N24 = 0; int N25 = 0;
	int N31 = 0; int N32 = 0; int N33 = 0; int N34 = 0; int N35 = 0;
	int N41 = 0; int N42 = 0; int N43 = 0; int N44 = 0; int N45 = 0;
	int N51 = 0; int N52 = 0; int N53 = 0; int N54 = 0; int N55 = 0;

	for (int i = 0; i < ref_Xsize * ref_Ysize; i++)
	{
		if ((int)_refBand[i] == 0)
		{
			if (_resBand[i] == 0)
			{
				N11 += 1;
			}
			else if (_resBand[i] == 50)
			{
				N12 += 1;
			}
			else if (_resBand[i] == 100)
			{
				N13 += 1;
			}
			else if (_resBand[i] == 150)
			{
				N14 += 1;
			}
			else
			{
				N15 += 1;
			}
		}
		else if ((int)_refBand[i] == 50)
		{
			if (_resBand[i] == 0)
			{
				N21 += 1;
			}
			else if (_resBand[i] == 50)
			{
				N22 += 1;
			}
			else if (_resBand[i] == 100)
			{
				N23 += 1;
			}
			else if (_resBand[i] == 150)
			{
				N24 += 1;
			}
			else
			{
				N25 += 1;
			}
		}
		else if ((int)_refBand[i] == 100)
		{
			if (_resBand[i] == 0)
			{
				N31 += 1;
			}
			else if (_resBand[i] == 50)
			{
				N32 += 1;
			}
			else if (_resBand[i] == 100)
			{
				N33 += 1;
			}
			else if (_resBand[i] == 150)
			{
				N34 += 1;
			}
			else
			{
				N35 += 1;
			}
		}
		else if ((int)_refBand[i] == 150)
		{
			if (_resBand[i] == 0)
			{
				N41 += 1;
			}
			else if (_resBand[i] == 50)
			{
				N42 += 1;
			}
			else if (_resBand[i] == 100)
			{
				N43 += 1;
			}
			else if (_resBand[i] == 150)
			{
				N44 += 1;
			}
			else
			{
				N45 += 1;
			}
		}
		else
		{
			if (_resBand[i] == 0)
			{
				N51 += 1;
			}
			else if (_resBand[i] == 50)
			{
				N52 += 1;
			}
			else if (_resBand[i] == 100)
			{
				N53 += 1;
			}
			else if (_resBand[i] == 150)
			{
				N54 += 1;
			}
			else
			{
				N55 += 1;
			}
		}

	}

	// ������ȷ��
	double c1 = (double)N11 / (N11 + N12 + N13 + N14 + N15);
	double c2 = (double)N22 / (N21 + N22 + N23 + N24 + N25);
	double c3 = (double)N33 / (N31 + N32 + N33 + N34 + N35);
	double c4 = (double)N44 / (N41 + N42 + N43 + N44 + N45);
	double c5 = (double)N55 / (N51 + N52 + N53 + N54 + N55);
	double c = (c1 + c2 + c3 + c4 + c5) / 5;
	cout << "��1����ȷ��Ϊ��" << c1 << endl;
	cout << "��2����ȷ��Ϊ��" << c2 << endl;
	cout << "��3����ȷ��Ϊ��" << c3 << endl;
	cout << "��4����ȷ��Ϊ��" << c4 << endl;
	cout << "��5����ȷ��Ϊ��" << c5 << endl;
	cout << "��ȷ�ʾ�ֵ��" << c << endl;


	//��������
	//int A_refChange = N11 + N21;
	//int B_refNoChange = N12 + N22;
	//int C_resChange = N11 + N12;
	//int D_resNoChange = N21 + N22;
	//int n1 = A_refChange + B_refNoChange;
	//int n2 = C_resChange + D_resNoChange;
	//int T_Count = N11 + N12 + N21 + N22;

	//����©���ʡ��龯�ʡ��ܷ��ྫ�ȡ�Kappaϵ��
	//double p1 = (double)N21 / A_refChange;
	//double p2 = (double)N12 / C_resChange;
	//double p = (double)(N11 + N22) / T_Count;
	//double kappa = (double)(T_Count * (N11 + N22) - (A_refChange * C_resChange + B_refNoChange * D_resNoChange)) /
	//	(T_Count * T_Count - (A_refChange * C_resChange + B_refNoChange * D_resNoChange));

	std::cout << "����������ɣ����������dataĿ¼��" << std::endl;

	std::fstream f;
	f.open(".\\data\\EvaResult.txt", std::ios::out);
	//д�������
	//f << "©���ʣ�" << p1 * 100 << "%" << '\n'
	//	<< "�龯�ʣ�" << p2 * 100 << "%" << '\n'
	//	<< "�ܷ��ྫ�ȣ�" << p * 100 << "%" << '\n'
	//	<< "Kappaϵ����" << kappa << '\n';

}

// ���ͼ��
void saveResult(Data& Img, int* predictLabels, const char* resultPath)
{

	//��ȡ����
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//���岨������˳��
	int Bandmap[1] = { 1 };

	int width = Img.getXsize();
	int height = Img.getYsize();

	//��������Ӱ������ݼ�
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* imgDstDS = Driver->Create(resultPath, width, height, 1, GDT_Int32, papszOptions);
	//������ʧ��
	if (!imgDstDS)
	{
		cout << "Create file failed!" << endl;
	}
	else
	{
		//д��ͼ��
		imgDstDS->RasterIO(GF_Write, 0, 0, width, height, predictLabels, width, height, GDT_Int32, 1, Bandmap, 0, 0, 0);
		cout << "����ɹ�����鿴�ļ���" << endl;
	}

	GDALClose(imgDstDS);
}

// �����Ȼ������
int maxLikelihoodClassify(vector<double> probbly, vector<Pixel> means, vector<Pixel> Var, int trainLabels[], int numTrain, Pixel testPix)
{
	int maxId = 0;
	double maxProb = 0;


	for (int i = 0; i < numTrain; i++)
	{

		double mean_r = means[i].r;
		double mean_g = means[i].g;
		double mean_b = means[i].b;
		double var_r = Var[i].r;
		double var_g = Var[i].g;
		double var_b = Var[i].b;

		double prob_r = gaussProDen(mean_r, var_r, testPix.r);
		double prob_g = gaussProDen(mean_g, var_g, testPix.g);
		double prob_b = gaussProDen(mean_b, var_b, testPix.b);
		double prob = prob_r * prob_g * prob_b * probbly[i];
		if (prob > maxProb)
		{
			maxProb = prob;
			maxId = trainLabels[i];
		}

	}
	return maxId;
}