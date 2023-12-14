#pragma once
#include "pch.h"
#include <iostream>
#include <gdal_priv.h>
#include <Eigen/Dense>
#include <fstream>
using namespace std;
using namespace Eigen;

class PointData
{
public:
	PointData(const char*);//有参构造，传入影像的地址进行初始化
	~PointData();
	int getXsize();
	int getBandnum();
	int getYsize();
	

	GDALDataset* Dataset;
	GDALDataType DataType;
	int Xsize;
	int Ysize;
	int Bandnum;

};

//计算一次多项式
tuple<Vector3d, Vector3d, Vector3d, Vector3d> PolynomialCal(const char*);

//几何校正
void Geometric_Correction(PointData&, PointData&, const char*, const char*);

//NCC_match
void Match_Ncc(PointData&, PointData&);