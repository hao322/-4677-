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
	PointData(const char*);//�вι��죬����Ӱ��ĵ�ַ���г�ʼ��
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

//����һ�ζ���ʽ
tuple<Vector3d, Vector3d, Vector3d, Vector3d> PolynomialCal(const char*);

//����У��
void Geometric_Correction(PointData&, PointData&, const char*, const char*);

//NCC_match
void Match_Ncc(PointData&, PointData&);