#include "pch.h"
#include "PointData.h"


int main()
{
	//注册文件格式
	GDALAllRegister();
	//设置支持中文路径
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* RefImg = ".\\data\\参考影像.tif";
	const char* CorImg = ".\\data\\待配准.tif";

	PointData* Cordata = new PointData(CorImg);
	PointData* Refdata = new PointData(RefImg);

	const char* pdata = ".\\data\\CTP.txt";
	const char* resultPath = ".\\data\\result.tif";

	//Geometric_Correction(*Refdata, *Cordata, pdata, resultPath);
	//Match_Ncc(*Cordata, *Refdata);

	system("pause");
	return 0;
}

