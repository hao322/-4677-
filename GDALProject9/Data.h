#pragma once
#include <gdal_priv.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

struct Pixel
{
	double r, g, b;
};
class Data
{
public:
	Data(const char*);
	int getXsize();
	int getYsize();
	int getBandnum();
	GDALDataset* getDataset();
	GDALDataType getDatatype();
	~Data();

private:
	GDALDataset* Dataset;
	GDALDataType Datatype;
	int Xsize;
	int Ysize;
	int Bandnum; 
};

double caldistance(Pixel, Pixel);
// 加载图像信息
vector<Pixel> loadImage(Data&);
// 最小距离分类
int minDistClassify(vector<Pixel>, int[], int, Pixel);
// 获取概率
vector<double> getProbability(vector<Pixel>, vector<Pixel>);
// 最大似然分类
int maxLikelihoodClassify(vector<double>, vector<Pixel>, vector<Pixel>, int[], int, Pixel);
// 保存影像
void saveResult(Data&, int*, const char*);

// 求均值及协方差
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeans(vector<Pixel>, vector<Pixel>);

// 精度评价
void EvaAccuracy(Data&, Data&);

// 像素2向量
Vector3d Pix2Vec(Pixel);

// 多元高斯概率密度计算
double testgauss(VectorXd, MatrixXd, VectorXd);

// 最大似然法测试
int testMax(vector<double>, vector<Vector3d>, vector<Matrix3d>, int trainLabels[], int numTrain, Pixel testPix);




