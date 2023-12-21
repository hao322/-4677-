#pragma once
#include <gdal_priv.h>
#include <iostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <Eigen/Dense>
#include <iomanip>
using namespace std;
using namespace Eigen;

// 图像的宽度和高度
const int Width = 500;
const int Height = 1000;

const double PI = 4 * atan(1);

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


// 加载图像信息
vector<Pixel> loadImage(Data&);
// 计算向量距离
double caldistance(Vector3d, Vector3d);
// 像素To向量
Vector3d Pix2Vec(Pixel);
// 获取均值和协方差矩阵
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeansCovPro(vector<Pixel>, vector<Pixel>);
// 多元高斯概率密度计算
double MulGaussProbability(VectorXd, MatrixXd, VectorXd);

// 最小距离法分类
int minDistClassify(vector<Vector3d>, int[], int, Pixel);
// 最大似然法分类
int maxLikelihoodClassify(vector<double>, vector<Vector3d>, vector<Matrix3d>, int trainLabels[], int numTrain, Pixel testPix);
// 保存影像
void saveResult(Data&, int*, const char*);

// 混淆矩阵计算accuracy及kappa系数
tuple<double, double> getAccuracyKappa(MatrixXi);
// 精度评价
void EvaAccuracy(vector<Pixel>, int[]);


