#pragma once
#include <gdal_priv.h>
#include <iostream>
#include <vector>
#include <random>
using namespace std;

class DataIni
{
public:
	DataIni(const char*);
	int getXsize();
	int getYsize();
	int getBandnum();
	GDALDataset* getDataset();
	GDALDataType getDatatype();
	~DataIni();

private:
	GDALDataset* Dataset;
	GDALDataType Datatype;
	int Xsize;
	int Ysize;
	int Bandnum;
};

