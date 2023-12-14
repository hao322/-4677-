#include "pch.h"
#include "DataIni.h"

DataIni::DataIni(const char* FilePath)
{
    this->Dataset = (GDALDataset*)GDALOpen(FilePath, GA_ReadOnly);
    this->Datatype = Dataset->GetRasterBand(1)->GetRasterDataType();
    this->Xsize = Dataset->GetRasterXSize();
    this->Ysize = Dataset->GetRasterYSize();
    this->Bandnum = Dataset->GetRasterCount();
}
int DataIni::getXsize()
{
    return Xsize;
}
int DataIni::getYsize()
{
    return Ysize;
}
int DataIni::getBandnum()
{
    return Bandnum;
}
GDALDataset* DataIni::getDataset()
{
    return Dataset;
}
GDALDataType DataIni::getDatatype()
{
    return Datatype;
}
DataIni::~DataIni() {};
