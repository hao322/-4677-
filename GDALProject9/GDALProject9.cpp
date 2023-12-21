#include "pch.h"
#include "Data.h"


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
	// 读取测试参考图像
	vector<Pixel> teLabelPl = loadImage(*testLabel);

	// 定义接受变量
	// 均值
	vector<Vector3d> means;
	// 协方差矩阵
	vector<Matrix3d> covMat;
	// 先验概率
	vector<double> probbly;
	tie(means, covMat, probbly) = getMeansCovPro(trPl, trLabelPl);

	// 预测每个像素
	int* minpredictLabels = new int[Width * Height];
	int* maxpredictLabels = new int[Width * Height];

	// 逐像素遍历
	for (int i = 0; i < Width * Height; i++)
	{
		//minpredictLabels[i] = minDistClassify(means, trainLabels, numTrain, testPl[i]);
		maxpredictLabels[i] = maxLikelihoodClassify(probbly, means, covMat, trainLabels, numTrain, testPl[i]);
	}

	const char* minresult = ".\\data\\最小距离法result_4677.tiff";
	const char* maxresult = ".\\data\\最大似然法result_4677.tiff";
	// 保存结果图像
	
	// 最小距离法
	//saveResult(*testImg, minpredictLabels, minresult);
	// 最大似然法
	//saveResult(*trainImg, maxpredictLabels, maxresult);


	//EvaAccuracy(teLabelPl, minpredictLabels);
	EvaAccuracy(teLabelPl, maxpredictLabels);

    system("pause");
    return -1;
}

