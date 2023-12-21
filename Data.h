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

// ͼ��Ŀ�Ⱥ͸߶�
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


// ����ͼ����Ϣ
vector<Pixel> loadImage(Data&);
// ������������
double caldistance(Vector3d, Vector3d);
// ����To����
Vector3d Pix2Vec(Pixel);
// ��ȡ��ֵ��Э�������
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeansCovPro(vector<Pixel>, vector<Pixel>);
// ��Ԫ��˹�����ܶȼ���
double MulGaussProbability(VectorXd, MatrixXd, VectorXd);

// ��С���뷨����
int minDistClassify(vector<Vector3d>, int[], int, Pixel);
// �����Ȼ������
int maxLikelihoodClassify(vector<double>, vector<Vector3d>, vector<Matrix3d>, int trainLabels[], int numTrain, Pixel testPix);
// ����Ӱ��
void saveResult(Data&, int*, const char*);

// �����������accuracy��kappaϵ��
tuple<double, double> getAccuracyKappa(MatrixXi);
// ��������
void EvaAccuracy(vector<Pixel>, int[]);


