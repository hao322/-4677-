#include "pch.h"
#include "Data.h"

// 图像的宽度和高度
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

// 计算两个向量距离
double caldistance(Pixel a, Pixel b)
{
	return sqrt(pow(a.r - b.r, 2) + pow(a.g - b.g, 2) + pow(a.b - b.b, 2));
}
// 读取图像数据
vector<Pixel> loadImage(Data& Img)
{
	using namespace std;
	vector<Pixel> pl;
	GDALDataType ImgType = Img.getDatatype();
	// 样本数据
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

// 样本均值1、样本均值2、协方差、先验概率
tuple<vector<Pixel>, vector<Vector3d>,vector<Matrix3d>, vector<double> > getMeans(vector<Pixel> train, vector<Pixel> Label)
{
	// 均值1、以结构体存储均值，用于最小距离法
	vector<Pixel> meanPixels(5);
	int k0 = 0; int k1 = 0; int k2 = 0; int k3 = 0; int k4 = 0;
	MatrixXd ccvv0(Height*Width, 3);
	MatrixXd ccvv1(Height*Width, 3);
	MatrixXd ccvv2(Height*Width, 3);
	MatrixXd ccvv3(Height*Width, 3);
	MatrixXd ccvv4(Height*Width, 3);
	for (int i = 0; i < Width * Height; i++)
	{
		if (Label[i].r == 0)
		{
			ccvv0(k0, 0) = train[i].r;
			ccvv0(k0, 1) = train[i].g;
			ccvv0(k0, 2) = train[i].b;
			k0++;
		}
		else if (Label[i].r == 50)
		{
			ccvv1(k1, 0) = train[i].r;
			ccvv1(k1, 1) = train[i].g;
			ccvv1(k1, 2) = train[i].b;
			k1++;
		}
		else if (Label[i].r == 100)
		{
			ccvv2(k2, 0) = train[i].r;
			ccvv2(k2, 1) = train[i].g;
			ccvv2(k2, 2) = train[i].b;
			k2++;
		}
		else if (Label[i].r == 150)
		{
			ccvv3(k3, 0) = train[i].r;
			ccvv3(k3, 1) = train[i].g;
			ccvv3(k3, 2) = train[i].b;
			k3++;
		}
		else
		{
			ccvv4(k4, 0) = train[i].r;
			ccvv4(k4, 1) = train[i].g;
			ccvv4(k4, 2) = train[i].b;
			k4++;
		}
	}
	ccvv0.resize(k0, 3); ccvv1.resize(k1, 3); ccvv2.resize(k2, 3); ccvv3.resize(k3, 3); ccvv4.resize(k4, 3);

	

	//double sumR0 = 0;	double sumG0 = 0;	double sumB0 = 0;
	//double sumR50 = 0;	double sumG50 = 0;	double sumB50 = 0;
	//double sumR100 = 0;	double sumG100 = 0;	double sumB100 = 0;
	//double sumR150 = 0;	double sumG150 = 0;	double sumB150 = 0;
	//double sumR200 = 0;	double sumG200 = 0;	double sumB200 = 0;

	//for (int mark = 0; mark < Height*Width; mark++)
	//{
	//	if (Label[mark].r == 0)
	//	{
	//		sumR0 += train[mark].r;
	//		sumG0 += train[mark].g;
	//		sumB0 += train[mark].b;
	//		k0++;
	//	}
	//	else if (Label[mark].r == 50)
	//	{
	//		sumR50 += train[mark].r;
	//		sumG50 += train[mark].g;
	//		sumB50 += train[mark].b;
	//		k1++;
	//	}
	//	else if (Label[mark].r == 100)
	//	{
	//		sumR100 += train[mark].r;
	//		sumG100 += train[mark].g;
	//		sumB100 += train[mark].b;
	//		k2++;
	//	}
	//	else if (Label[mark].r == 150)
	//	{
	//		sumR150 += train[mark].r;
	//		sumG150 += train[mark].g;
	//		sumB150 += train[mark].b;
	//		k3++;
	//	}
	//	else
	//	{
	//		sumR200 += train[mark].r;
	//		sumG200 += train[mark].g;
	//		sumB200 += train[mark].b;
	//		k4++;
	//	}
	//}

	// 先验概率
	vector<double> probability(5);
	probability[0] = (double)k0 / train.size();
	probability[1] = (double)k1 / train.size();
	probability[2] = (double)k2 / train.size();
	probability[3] = (double)k3 / train.size();
	probability[4] = (double)k4 / train.size();

	//meanPixels[0].r = sumR0 / k0; meanPixels[0].g = sumG0 / k0;	meanPixels[0].b = sumB0 / k0;
	//meanPixels[1].r = sumR50 / k1; meanPixels[1].g = sumG50 / k1; meanPixels[1].b = sumB50 / k1;
	//meanPixels[2].r = sumR100 / k2; meanPixels[2].g = sumG100 / k2;	meanPixels[2].b = sumB100 / k2;
	//meanPixels[3].r = sumR150 / k3; meanPixels[3].g = sumG150 / k3;	meanPixels[3].b = sumB150 / k3;
	//meanPixels[4].r = sumR200 / k4;  meanPixels[4].g = sumG200 / k4; meanPixels[4].b = sumB200 / k4;

	// 均值2、以向量存储均值，用于最大似然法
	vector<Vector3d> means;
	MatrixXd a0 = ccvv0.colwise().mean();
	a0.resize(a0.rows(), a0.cols());
	Vector3d m0 = a0.transpose();
	MatrixXd a1 = ccvv1.colwise().mean();
	a1.resize(a1.rows(), a1.cols());
	Vector3d m1 = a1.transpose();
	MatrixXd a2 = ccvv2.colwise().mean();
	a2.resize(a2.rows(), a2.cols());
	Vector3d m2 = a2.transpose();
	MatrixXd a3 = ccvv3.colwise().mean();
	a3.resize(a3.rows(), a3.cols());
	Vector3d m3 = a3.transpose();
	MatrixXd a4 = ccvv4.colwise().mean();
	a4.resize(a4.rows(), a4.cols());
	Vector3d m4 = a4.transpose();

	means.push_back(m0);
	means.push_back(m1);
	means.push_back(m2);
	means.push_back(m3);
	means.push_back(m4);

	//means[0] = ccvv0.colwise().mean().transpose();
	//means[1] = ccvv1.colwise().mean().transpose();
	//means[2] = ccvv2.colwise().mean().transpose();
	//means[3] = ccvv3.colwise().mean().transpose();
	//means[4] = ccvv4.colwise().mean().transpose();

	MatrixXd cov_0(ccvv0.rows(), 3);
	MatrixXd cov_1(ccvv1.rows(), 3);
	MatrixXd cov_2(ccvv2.rows(), 3);
	MatrixXd cov_3(ccvv3.rows(), 3);
	MatrixXd cov_4(ccvv4.rows(), 3);

	//1 求各个维度均值
	//MatrixXd mean_0 = ccvv0.colwise().mean();
	//MatrixXd mean_1 = ccvv0.colwise().mean();
	//MatrixXd mean_2 = ccvv0.colwise().mean();
	//MatrixXd mean_3 = ccvv0.colwise().mean();

	//2 每个样本减去均值
	for (int i = 0; i < ccvv0.rows(); i++)
	{
		cov_0.row(i) = ccvv0.row(i) - a0;
	}
	for (int i = 0; i < ccvv1.rows(); i++)
	{
		cov_1.row(i) = ccvv1.row(i) - a1;
	}
	for (int i = 0; i < ccvv2.rows(); i++)
	{
		cov_2.row(i) = ccvv2.row(i) - a2;
	}
	for (int i = 0; i < ccvv3.rows(); i++)
	{
		cov_3.row(i) = ccvv3.row(i) - a3;
	}
	for (int i = 0; i < ccvv4.rows(); i++)
	{
		cov_4.row(i) = ccvv4.row(i) - a4;
	}
	//3 计算
	cov_0 = cov_0.transpose() * cov_0 / (cov_0.rows() - 1);
	cov_1 = cov_1.transpose() * cov_1 / (cov_1.rows() - 1);
	cov_2 = cov_2.transpose() * cov_2 / (cov_2.rows() - 1);
	cov_3 = cov_3.transpose() * cov_3 / (cov_3.rows() - 1);
	cov_4 = cov_4.transpose() * cov_4 / (cov_4.rows() - 1);


	//for (int mark = 0; mark < Height*Width; mark++)
	//{
	//	if (Label[mark].r == 0)
	//	{
	//		cov0(a0, 0) = train[mark].r - meanPixels[0].r;
	//		cov0(a0, 1) = train[mark].g - meanPixels[0].g;
	//		cov0(a0, 2) = train[mark].b - meanPixels[0].b;
	//		//cov0(a0, 0) = train[mark].r;
	//		//cov0(a0, 1) = train[mark].g;
	//		//cov0(a0, 2) = train[mark].b;
	//		a0++;
	//	}
	//	else if (Label[mark].r == 50)
	//	{
	//		cov1(a1, 0) = train[mark].r - meanPixels[1].r;
	//		cov1(a1, 1) = train[mark].g - meanPixels[1].g;
	//		cov1(a1, 2) = train[mark].b - meanPixels[1].b;
	//		//cov1(a1, 0) = train[mark].r;
	//		//cov1(a1, 1) = train[mark].g;
	//		//cov1(a1, 2) = train[mark].b;
	//		a1++;
	//	}
	//	else if (Label[mark].r == 100)
	//	{
	//		cov2(a2, 0) = train[mark].r - meanPixels[2].r;
	//		cov2(a2, 1) = train[mark].g - meanPixels[2].g;
	//		cov2(a2, 2) = train[mark].b - meanPixels[2].b;
	//		//cov2(a2, 0) = train[mark].r;
	//		//cov3(a2, 1) = train[mark].g;
	//		//cov2(a2, 2) = train[mark].b;
	//		a2++;
	//	}
	//	else if (Label[mark].r == 150)
	//	{
	//		cov3(a3, 0) = train[mark].r - meanPixels[3].r;
	//		cov3(a3, 1) = train[mark].g - meanPixels[3].g;
	//		cov3(a3, 2) = train[mark].b - meanPixels[3].b;
	//		//cov3(a3, 0) = train[mark].r;
	//		//cov3(a3, 1) = train[mark].g;
	//		//cov3(a3, 2) = train[mark].b;
	//		a3++;
	//	}
	//	else
	//	{
	//		cov4(a4, 0) = train[mark].r - meanPixels[4].r;
	//		cov4(a4, 1) = train[mark].g - meanPixels[4].g;
	//		cov4(a4, 2) = train[mark].b - meanPixels[4].b;
	//		//cov4(a4, 0) = train[mark].r;
	//		//cov4(a4, 1) = train[mark].g;
	//		//cov4(a4, 2) = train[mark].b;
	//		a4++;
	//	}
	//}
	
	//Matrix3d Cov0 = cov0.transpose() * cov0 / (cov0.cols() - 1);
	//Matrix3d Cov1 = cov1.transpose() * cov1 / (cov1.cols() - 1);
	//Matrix3d Cov2 = cov2.transpose() * cov2 / (cov2.cols() - 1);
	//Matrix3d Cov3 = cov3.transpose() * cov3 / (cov3.cols() - 1);
	//Matrix3d Cov4 = cov4.transpose() * cov4 / (cov4.cols() - 1);
	
	vector<Matrix3d> Cov(5);
	Cov[0] = cov_0;
	Cov[1] = cov_1;
	Cov[2] = cov_2;
	Cov[3] = cov_3;
	Cov[4] = cov_4;
	//Cov[0] = getCovMat(cov0);
	//Cov[1] = getCovMat(cov1);
	//Cov[2] = getCovMat(cov2);
	//Cov[3] = getCovMat(cov3);
	//Cov[4] = getCovMat(cov4);

	// vector从0到4分别对应0、50、100、150、200
	return make_tuple(meanPixels, means, Cov, probability);
}

// 协方差矩阵计算
//Matrix3d getCovMat(MatrixXd samples)
//{
//	MatrixXd cov_;
//
//	//1 求各个维度均值
//	MatrixXd mean_ = samples.colwise().mean();
//
//	//2 每个样本减去均值
//	for (size_t i = 0; i < samples.rows(); i++)
//	{
//		cov_.row(i) = samples.row(i) - mean_;
//	}
//	//3 计算
//	cov_ = cov_.transpose() * cov_ / (cov_.rows() - 1);
//
//	return cov_;
//}

// 高斯概率密度计算
double gaussProDen(double means, double var, double x)
{
	double temp1 = -1 * (x - means) * (x - means) / (2 * var);
	double temp2 = 1 / sqrt(2 * PI * var);

	return temp2 * exp(temp1);
}

// 多元高斯概率密度计算
double testgauss(VectorXd mean, MatrixXd cov, VectorXd test)
{
	int k = mean.size();

	double determinant = cov.determinant();
	double normalization = 1.0 / (pow(2 * PI, k / 2) * sqrt(determinant));

	MatrixXd inv_cov = cov.inverse();
	double exponent = -0.5 * (test - mean).transpose() * inv_cov * (test - mean);
	return normalization * exp(exponent);
}

// 像素到向量的转换
Vector3d Pix2Vec(Pixel a)
{
	VectorXd V(3);
	V << a.r, a.g, a.b;
	return V;
}

// 最大似然法测试
int testMax(vector<double> probbly, vector<Vector3d> means, vector<Matrix3d> cov, int trainLabels[], int numTrain, Pixel testPix)
{
	Vector3d tePix = Pix2Vec(testPix);

	int maxId = -1;
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

// 最小距离分类
int minDistClassify(vector<Pixel> meanPix, int trainLabels[], int numTrain, Pixel testPix)
{
	double minDist = INT_MAX;
	unsigned int predLabel;

	for (int i = 0; i < numTrain; ++i)
	{
		// 这里应该是与各个均值中心之间的距离
		double dist = caldistance(testPix, meanPix[i]);

		if (dist < minDist)
		{
			minDist = dist;
			predLabel = trainLabels[i];
		}
	}
	// 预测类别
	return predLabel;
}

//精度评定
void EvaAccuracy(Data& ref, Data& res)
{
	//reference
	GDALDataType ref_type = ref.getDatatype(); //GDT_UInt16 对应 unsigned short
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

	// 计算正确率
	double c1 = (double)N11 / (N11 + N12 + N13 + N14 + N15);
	double c2 = (double)N22 / (N21 + N22 + N23 + N24 + N25);
	double c3 = (double)N33 / (N31 + N32 + N33 + N34 + N35);
	double c4 = (double)N44 / (N41 + N42 + N43 + N44 + N45);
	double c5 = (double)N55 / (N51 + N52 + N53 + N54 + N55);
	double c = (c1 + c2 + c3 + c4 + c5) / 5;
	cout << "第1类正确率为：" << c1 << endl;
	cout << "第2类正确率为：" << c2 << endl;
	cout << "第3类正确率为：" << c3 << endl;
	cout << "第4类正确率为：" << c4 << endl;
	cout << "第5类正确率为：" << c5 << endl;
	cout << "正确率均值：" << c << endl;


	//混淆矩阵
	//int A_refChange = N11 + N21;
	//int B_refNoChange = N12 + N22;
	//int C_resChange = N11 + N12;
	//int D_resNoChange = N21 + N22;
	//int n1 = A_refChange + B_refNoChange;
	//int n2 = C_resChange + D_resNoChange;
	//int T_Count = N11 + N12 + N21 + N22;

	//计算漏检率、虚警率、总分类精度、Kappa系数
	//double p1 = (double)N21 / A_refChange;
	//double p2 = (double)N12 / C_resChange;
	//double p = (double)(N11 + N22) / T_Count;
	//double kappa = (double)(T_Count * (N11 + N22) - (A_refChange * C_resChange + B_refNoChange * D_resNoChange)) /
	//	(T_Count * T_Count - (A_refChange * C_resChange + B_refNoChange * D_resNoChange));

	std::cout << "精度评定完成，结果保存于data目录下" << std::endl;

	std::fstream f;
	f.open(".\\data\\EvaResult.txt", std::ios::out);
	//写入的内容
	//f << "漏检率：" << p1 * 100 << "%" << '\n'
	//	<< "虚警率：" << p2 * 100 << "%" << '\n'
	//	<< "总分类精度：" << p * 100 << "%" << '\n'
	//	<< "Kappa系数：" << kappa << '\n';

}

// 输出图像
void saveResult(Data& Img, int* predictLabels, const char* resultPath)
{

	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//定义波段排列顺序
	int Bandmap[1] = { 1 };

	int width = Img.getXsize();
	int height = Img.getYsize();

	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* imgDstDS = Driver->Create(resultPath, width, height, 1, GDT_Int32, papszOptions);
	//若创建失败
	if (!imgDstDS)
	{
		cout << "Create file failed!" << endl;
	}
	else
	{
		//写入图像
		imgDstDS->RasterIO(GF_Write, 0, 0, width, height, predictLabels, width, height, GDT_Int32, 1, Bandmap, 0, 0, 0);
		cout << "保存成功，请查看文件！" << endl;
	}

	GDALClose(imgDstDS);
}

// 最大似然法分类
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
