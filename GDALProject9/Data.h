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
// ����ͼ����Ϣ
vector<Pixel> loadImage(Data&);
// ��С�������
int minDistClassify(vector<Pixel>, int[], int, Pixel);
// ��ȡ����
vector<double> getProbability(vector<Pixel>, vector<Pixel>);
// �����Ȼ����
int maxLikelihoodClassify(vector<double>, vector<Pixel>, vector<Pixel>, int[], int, Pixel);
// ����Ӱ��
void saveResult(Data&, int*, const char*);

// ���ֵ��Э����
tuple<vector<Vector3d>, vector<Matrix3d>, vector<double> > getMeans(vector<Pixel>, vector<Pixel>);

// ��������
void EvaAccuracy(Data&, Data&);

// ����2����
Vector3d Pix2Vec(Pixel);

// ��Ԫ��˹�����ܶȼ���
double testgauss(VectorXd, MatrixXd, VectorXd);

// �����Ȼ������
int testMax(vector<double>, vector<Vector3d>, vector<Matrix3d>, int trainLabels[], int numTrain, Pixel testPix);




