#include "pch.h"
#include "Data.h"

const int Xsize = 500;
const int Ysize = 1000;

// 5个类别
const int numTrain = 5;
int trainLabels[5] = { 0, 50, 100, 150, 200 };

int main()
{
	//驱动设置和中文路径设置
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

	const char* trImg = ".\\data\\train.tiff";
	const char* trLabel = ".\\data\\train_label.tiff";
	const char* teImg = ".\\data\\test.tiff";
	const char* teLabel = ".\\data\\test_label.tiff";

	Data* trainImg = new Data(trImg);
	Data* trainLabel = new Data(trLabel);
	Data* testImg = new Data(teImg);
	Data* testLabel = new Data(teLabel);
	
	// 读取训练图像
	vector<Pixel> trPl = loadImage(*trainImg);
	// 读取训练标签图像
	vector<Pixel> trLabelPl = loadImage(*trainLabel);
	// 读取测试图像
	vector<Pixel> testPl = loadImage(*testImg);


	//	// 利用tie进行解包元素的值
	// 定义接受变量
	// 均值
	vector<Vector3d> meansVec;
	// 协方差
	vector<Matrix3d> var;
	// 先验概率
	vector<double> probbly;
	tie(meansVec, var, probbly) = getMeans(trPl, trLabelPl);

	// 预测每个像素
	int* minpredictLabels = new int[Xsize * Ysize];
	int* maxpredictLabels = new int[Xsize * Ysize];
	int flag = 0;
	// 逐像素遍历
	for (int i = 0; i < Ysize; i++)
	{
		for (int j = 0; j < Xsize; j++)
		{
			flag = i * Xsize + j;
			Pixel tePix = testPl[flag];
			// minpredictLabels[flag] = minDistClassify(means, trainLabels, numTrain, tePix);
			// maxpredictLabels[flag] = maxLikelihoodClassify(probbly, meanPix, varPix, trainLabels, numTrain, tePix);
			maxpredictLabels[flag] = testMax(probbly, meansVec, var, trainLabels, numTrain, tePix);
		}
	}

	const char* minresult = ".\\data\\result_最小距离法.tiff";
	const char* maxresult = ".\\data\\result_最大似然法_测试新.tiff";
	// 保存结果图像
	// saveResult(*trainImg, minpredictLabels, minresult);
	saveResult(*trainImg, maxpredictLabels, maxresult);

	//Data* mintest = new Data(minresult);
	//Data* maxtest = new Data(maxresult);
	//cout << "最小距离法" << endl;
	//EvaAccuracy(*testLabel, *mintest);
	//cout << "最大似然法" << endl;
	//EvaAccuracy(*testLabel, *maxtest);


    system("pause");
    return -1;
}

