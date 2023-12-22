#include "pch.h"
#include "Data.h"

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
double caldistance(Vector3d a, Vector3d b)
{
	return sqrt(pow(a[0] - b[0], 2) + pow(a[1] - b[1], 2) + pow(a[2] - b[2], 2));
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
// 像素到向量的转换
Vector3d Pix2Vec(Pixel a)
{
	VectorXd V(3);
	V << a.r, a.g, a.b;
	return V;
}

//获取均值、协方差矩阵、先验概率
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeansCovPro(vector<Pixel> train, vector<Pixel> Label)
{
	int k0 = 0; int k1 = 0; int k2 = 0; int k3 = 0; int k4 = 0;
	MatrixXd ccvv0(Height * Width, 3);
	MatrixXd ccvv1(Height * Width, 3);
	MatrixXd ccvv2(Height * Width, 3);
	MatrixXd ccvv3(Height * Width, 3);
	MatrixXd ccvv4(Height * Width, 3);
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

	ccvv0.conservativeResize(k0, 3); 
	ccvv1.conservativeResize(k1, 3); 
	ccvv2.conservativeResize(k2, 3); 
	ccvv3.conservativeResize(k3, 3); 
	ccvv4.conservativeResize(k4, 3);

	// 先验概率
	vector<double> probability(5);
	probability[0] = (double)k0 / train.size();
	probability[1] = (double)k1 / train.size();
	probability[2] = (double)k2 / train.size();
	probability[3] = (double)k3 / train.size();
	probability[4] = (double)k4 / train.size();

	// 以向量存储均值，用于最大似然法
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

	MatrixXd cov_0(ccvv0.rows(), 3);
	MatrixXd cov_1(ccvv1.rows(), 3);
	MatrixXd cov_2(ccvv2.rows(), 3);
	MatrixXd cov_3(ccvv3.rows(), 3);
	MatrixXd cov_4(ccvv4.rows(), 3);

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

	vector<Matrix3d> Cov(5);
	Cov[0] = cov_0;
	Cov[1] = cov_1;
	Cov[2] = cov_2;
	Cov[3] = cov_3;
	Cov[4] = cov_4;

	// vector从0到4分别对应0、50、100、150、200
	return make_tuple(means, Cov, probability);


}

// 多元高斯概率密度计算
double MulGaussProbability(VectorXd mean, MatrixXd cov, VectorXd test)
{
	int k = mean.size();

	double determinant = cov.determinant();
	double normalization = 1.0 / (pow(2 * PI, k / 2) * sqrt(determinant));

	MatrixXd inv_cov = cov.inverse();
	double exponent = -0.5 * (test - mean).transpose() * inv_cov * (test - mean);
	return normalization * exp(exponent);
}

// 最大似然法分类
int maxLikelihoodClassify(vector<double> probbly, vector<Vector3d> means, vector<Matrix3d> cov, int trainLabels[], int numTrain, Pixel testPix)
{
	Vector3d tePix = Pix2Vec(testPix);

	int maxId = -1;
	double maxProb = 0;

	for (int i = 0; i < numTrain; i++)
	{
		double prob = MulGaussProbability(means[i], cov[i], tePix) * probbly[i];
		if (prob > maxProb)
		{
			maxProb = prob;
			maxId = trainLabels[i];
		}
	}
	return maxId;
}

// 最小距离分类
int minDistClassify(vector<Vector3d> meanPix, int trainLabels[], int numTrain, Pixel testPix)
{
	Vector3d testVec = Pix2Vec(testPix);
	double minDist = INT_MAX;
	unsigned int predLabel;

	for (int i = 0; i < numTrain; ++i)
	{
		// 这里应该是与各个均值中心之间的距离
		double dist = caldistance(testVec, meanPix[i]);

		if (dist < minDist)
		{
			minDist = dist;
			predLabel = trainLabels[i];
		}
	}
	// 预测类别
	return predLabel;
}

// 保存图像
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

//精度评定
void EvaAccuracy(vector<Pixel> reference, int* result)
{
	MatrixXi Nxx = MatrixXi::Zero(5, 5);

	// 计算各项
	for (int i = 0; i < Width * Height; i++)
	{
		if (reference[i].r==0)
		{
			if (result[i]==0)
			{
				Nxx(0, 0) += 1;
			}
			else if (result[i]==50)
			{
				Nxx(0, 1) += 1;
			}
			else if (result[i] == 100)
			{
				Nxx(0, 2) += 1;
			}
			else if (result[i] == 150)
			{
				Nxx(0, 3) += 1;
			}
			else
			{
				Nxx(0, 4) += 1;
			}

		}
		else if (reference[i].r==50)
		{
			if (result[i] == 0)
			{
				Nxx(1, 0) += 1;
			}
			else if (result[i] == 50)
			{
				Nxx(1, 1) += 1;
			}
			else if (result[i] == 100)
			{
				Nxx(1, 2) += 1;
			}
			else if (result[i] == 150)
			{
				Nxx(1, 3) += 1;
			}
			else
			{
				Nxx(1, 4) += 1;
			}
		}
		else if (reference[i].r == 100)
		{
			if (result[i] == 0)
			{
				Nxx(2, 0) += 1;
			}
			else if (result[i] == 50)
			{
				Nxx(2, 1) += 1;
			}
			else if (result[i] == 100)
			{
				Nxx(2, 2) += 1;
			}
			else if (result[i] == 150)
			{
				Nxx(2, 3) += 1;
			}
			else
			{
				Nxx(2, 4) += 1;
			}
		}
		else if (reference[i].r == 150)
		{
			if (result[i] == 0)
			{
				Nxx(3, 0) += 1;
			}
			else if (result[i] == 50)
			{
				Nxx(3, 1) += 1;
			}
			else if (result[i] == 100)
			{
				Nxx(3, 2) += 1;
			}
			else if (result[i] == 150)
			{
				Nxx(3, 3) += 1;
			}
			else
			{
				Nxx(3, 4) += 1;
			}
		}
		else
		{
			if (result[i] == 0)
			{
				Nxx(4, 0) += 1;
			}
			else if (result[i] == 50)
			{
				Nxx(4, 1) += 1;
			}
			else if (result[i] == 100)
			{
				Nxx(4, 2) += 1;
			}
			else if (result[i] == 150)
			{
				Nxx(4, 3) += 1;
			}
			else
			{
				Nxx(4, 4) += 1;
			}
		}
	}

	double sumCount = Nxx.sum();
	VectorXi rowCount = Nxx.rowwise().sum();
	VectorXi colCount = Nxx.colwise().sum();
	// 计算TP、FP、TN、FN
	int TP0 = Nxx(0, 0); int TP1 = Nxx(1, 1); int TP2 = Nxx(2, 2); int TP3 = Nxx(3, 3); int TP4 = Nxx(4, 4);

	//int FP0 = colCount[0] - Nxx(0, 0);
	//int FP1 = colCount[1] - Nxx(1, 1);
	//int FP2 = colCount[2] - Nxx(2, 2);
	//int FP3 = colCount[3] - Nxx(3, 3);
	//int FP4 = colCount[4] - Nxx(4, 4);

	// 计算精确率、召回率、总体精度
	double Precision_0 = (double)TP0 / colCount[0];
	double Recall_0 = (double)TP0 / rowCount[0];
	double Precision_1 = (double)TP1 / colCount[1];
	double Recall_1 = (double)TP1 / rowCount[1];
	double Precision_2 = (double)TP2 / colCount[2];
	double Recall_2 = (double)TP2 / rowCount[2];
	double Precision_3 = (double)TP3 / colCount[3];
	double Recall_3 = (double)TP3 / rowCount[3];
	double Precision_4 = (double)TP4 / colCount[4];
	double Recall_4 = (double)TP4 / rowCount[4];

	double accuracy = double(Nxx(0, 0) + Nxx(1, 1) + Nxx(2, 2) + Nxx(3, 3) + Nxx(4, 4)) / Nxx.sum();

	std::cout << "精度评定完成，结果保存于data目录下" << std::endl;

	std::fstream f;
	f.open(".\\data\\EvaResult.txt", std::ios::out);
	//写入的内容
	f << "混淆矩阵如下：" << '\n';
	for (int i = 0; i < Nxx.rows(); i++)
	{
		for (int j = 0; j < Nxx.cols(); j++)
		{
			f << setw(6) << Nxx(i, j)<<" ";
			if (j==Nxx.cols()-1)
			{
				f << '\n';
			}
		}
	}

	f << "第1类正确率为：" << Precision_0 << '\n'
		<< "第1类召回率为：" << Recall_0 << '\n'
		<< "第2类正确率为：" << Precision_1 << '\n'
		<< "第2类召回率为：" << Recall_1 << '\n'
		<< "第3类正确率为：" << Precision_2 << '\n'
		<< "第3类召回率为：" << Recall_2 << '\n'
		<< "第4类正确率为：" << Precision_3 << '\n'
		<< "第4类召回率为：" << Recall_3 << '\n'
		<< "第5类正确率为：" << Precision_4 << '\n'
		<< "第5类召回率为：" << Recall_4 << '\n'
		<< "总体精度Accuracy：" << accuracy << '\n';
}