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

// 样本均值 及各类别样本数量
tuple<vector<Pixel>, int, int, int, int, int> getMeans(vector<Pixel> train, vector<Pixel> Label)
{
	int mark = 0; int k0 = 0; int k50 = 0; int k100 = 0; int k150 = 0; int k200 = 0;

	double sumR0 = 0;	double sumG0 = 0;	double sumB0 = 0;
	double sumR50 = 0;	double sumG50 = 0;	double sumB50 = 0;
	double sumR100 = 0;	double sumG100 = 0;	double sumB100 = 0;
	double sumR150 = 0;	double sumG150 = 0;	double sumB150 = 0;
	double sumR200 = 0;	double sumG200 = 0;	double sumB200 = 0;

	double r0mean, r50mean, r100mean, r150mean, r200mean;
	double g0mean, g50mean, g100mean, g150mean, g200mean;
	double b0mean, b50mean, b100mean, b150mean, b200mean;

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
				k50++;
			}
			else if (Label[mark].r == 100)
			{
				sumR100 += train[mark].r;
				sumG100 += train[mark].g;
				sumB100 += train[mark].b;
				k100++;
			}
			else if (Label[mark].r == 150)
			{
				sumR150 += train[mark].r;
				sumG150 += train[mark].g;
				sumB150 += train[mark].b;
				k150++;
			}
			else
			{
				sumR200 += train[mark].r;
				sumG200 += train[mark].g;
				sumB200 += train[mark].b;
				k200++;
			}
		}
	}
	r0mean = sumR0 / k0; g0mean = sumG0 / k0; b0mean = sumB0 / k0;
	r50mean = sumR50 / k50; g50mean = sumG50 / k50; b50mean = sumB50 / k50;
	r100mean = sumR100 / k100; g100mean = sumG100 / k100; b100mean = sumB100 / k100;
	r150mean = sumR150 / k150; g150mean = sumG150 / k150; b150mean = sumB150 / k150;
	r200mean = sumR200 / k200; g200mean = sumG200 / k200; b200mean = sumB200 / k200;

	vector<Pixel> meanPixels(5);
	meanPixels[0].r = r0mean; meanPixels[0].g = g0mean;	meanPixels[0].b = b0mean;
	meanPixels[1].r = r50mean; meanPixels[1].g = g50mean;	meanPixels[1].b = b50mean;
	meanPixels[2].r = r100mean; meanPixels[2].g = g100mean;	meanPixels[2].b = b100mean;
	meanPixels[3].r = r150mean; meanPixels[3].g = g150mean;	meanPixels[3].b = b150mean;
	meanPixels[4].r = r200mean; meanPixels[4].g = g200mean;	meanPixels[4].b = b200mean;

	// meanPixels vector从0到4分别对应0、50、100、150、200
	return make_tuple(meanPixels, k0, k50, k100, k150, k200);
}

// 样本方差
vector<Pixel> Variance(vector<Pixel> train, vector<Pixel> Label)
{
	// 利用tie进行解包元素的值
	// 定义接受变量
	vector<Pixel> means; // 求方差需要用均值
	int a0, a50, a100, a150, a200; // 各样本数量
	tie(means, a0, a50, a100, a150, a200) = getMeans(train, Label);

	double sumR0 = 0;	double sumG0 = 0;	double sumB0 = 0;
	double sumR50 = 0;	double sumG50 = 0;	double sumB50 = 0;
	double sumR100 = 0;	double sumG100 = 0;	double sumB100 = 0;
	double sumR150 = 0;	double sumG150 = 0;	double sumB150 = 0;
	double sumR200 = 0;	double sumG200 = 0;	double sumB200 = 0;
	
	// 求方差
	int mark = 0;
	for (int i = 0; i < Height; i++)
	{
		for (int j = 0; j < Width; j++)
		{
			mark = i * Width + j;
			if (Label[mark].r == 0)
			{
				sumR0 += train[mark].r * train[mark].r;
				sumG0 += train[mark].g * train[mark].g;
				sumB0 += train[mark].b * train[mark].b;
			}
			else if (Label[mark].r == 50)
			{
				sumR50 += train[mark].r * train[mark].r;
				sumG50 += train[mark].g * train[mark].g;
				sumB50 += train[mark].b * train[mark].b;
			}
			else if (Label[mark].r == 100)
			{
				sumR100 += train[mark].r * train[mark].r;
				sumG100 += train[mark].g * train[mark].g;
				sumB100 += train[mark].b * train[mark].b;
			}
			else if (Label[mark].r == 150)
			{
				sumR150 += train[mark].r * train[mark].r;
				sumG150 += train[mark].g * train[mark].g;
				sumB150 += train[mark].b * train[mark].g;
			}
			else
			{
				sumR200 += train[mark].r * train[mark].r;
				sumG200 += train[mark].g * train[mark].g;
				sumB200 += train[mark].b * train[mark].b;
			}
		}
	}

	vector<Pixel> axx(5);
	axx[0].r = sumR0 / a0; axx[0].g = sumG0 / a0; axx[0].b = sumB0 / a0;
	axx[1].r = sumR50 / a50; axx[1].g = sumG50 / a50; axx[1].b = sumB50 / a50;
	axx[2].r = sumR100 / a100; axx[2].g = sumG100 / a100; axx[2].b = sumB100 / a100;
	axx[3].r = sumR150 / a150; axx[3].g = sumG150 / a150; axx[3].b = sumB150 / a150;
	axx[4].r = sumR200 / a200; axx[4].g = sumG200 / a200; axx[4].b = sumB200 / a200;

	vector<Pixel> Var;
	for (int i = 0; i < 5; i++)
	{
		Pixel pl;
		pl.r = axx[i].r - means[i].r * means[i].r;
		pl.g = axx[i].g - means[i].g * means[i].g;
		pl.b = axx[i].b - means[i].b * means[i].b;
		Var.push_back(pl);
	}

	return Var;
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

// 样本数量获取，求得概率 P_Yi
vector<double> getProbability(vector<Pixel> train, vector<Pixel> Label)
{
	int mark = 0; int num0 = 0; int num50 = 0; int num100 = 0; int num150 = 0; int num200 = 0;
	vector<double> probbly(5);

	for (int i = 0; i < Width * Height; i++)
	{
		if (Label[i].r == 0)
		{
			num0 += 1;
		}
		else if (Label[i].r == 50)
		{
			num50 += 1;
		}
		else if (Label[i].r == 100)
		{
			num100 += 1;
		}
		else if (Label[i].r == 150)
		{
			num150 += 1;
		}
		else
		{
			num200 += 1;
		}
	}
	mark = Width * Height;
	probbly[0] = (double)num0 / mark;probbly[1] = (double)num50 / mark; probbly[2] = (double)num100 / mark; probbly[3] = (double)num150 / mark; probbly[4] = (double)num200 / mark;

	return probbly;
}

// 高斯概率密度计算
double gaussProDen(double means, double var, double x)
{
	double temp = -1 * (x - means) * (x - means) / (2 * var);
	double result = 1 / (sqrt(2 * PI * var)) * exp(temp);
	return result;
}

// 最大似然法分类
int maxLikelihoodClassify(vector<double> probbly,vector<Pixel> means, vector<Pixel> Var, int trainLabels[], int numTrain, Pixel testPix)
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
		//double prob = sqrt(prob_r * prob_r + prob_g * prob_g + prob_b * prob_b) * probbly[i];
		double prob = prob_r * prob_g * prob_b * probbly[i];
		if (prob > maxProb)
		{
			maxProb = prob;
			maxId = trainLabels[i];
		}

	}
	return maxId;
}

//均值滤波
int* meanFilter(int *Input)
{
	//申请内存 一个波段
	int* Output = new int[Width * Height];

	for (int i = 0; i < Height; i++) {
		for (int j = 0; j < Width; j++) {
			if (i == 0 || j == 0 || i == Height - 1 || j == Width - 1) {
				Output[Width * i + j] = Input[Width * i + j];
			}
		}
	}

	int temp[9] = { 0 };//3*3的
	for (int i = 1; i < Height - 1; i++) {
		for (int j = 1; j < Width - 1; j++) {
			for (int k = 0; k < 9; k++) {
				*(temp + k) = Input[Width * (i + k / 3) + j + k % 3];
			}
			int sum = 0;
			for (int t = 0; t < 9; t++) {
				sum += *(temp + t);
			}

			Output[Width * i + j] = sum / 9;
		}
	}
	return Output;
}

// 中值滤波
int* medianFilter(int* Input)
{

	int* Output = new int[Width * Height];

	for (int i = 0; i < Height; i++) {
		for (int j = 0; j < Width; j++) {
			if (i == 0 || j == 0 || i == Height- 1 || j == Width - 1) {
				Output[Width * i + j] = Input[Width * i + j];
			}
		}
	}
	int temp[9] = { 0 };//3*3的
	for (int i = 1; i < Height - 1; i++) {
		for (int j = 1; j < Width - 1; j++) {
			for (int k = 0; k < 9; k++) {
				*(temp + k) = Input[Width * (i + (k / 3 - 1)) + j + k % 3 - 1];
			}
			sort(temp, temp + 9);
			Output[Width * i + j] = (int)*(temp + 4);
		}
	}

	return Output;
}

void saveResult(Data& Img, int *predictLabels, const char* resultPath) 
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


