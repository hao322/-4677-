#include "pch.h"
#include "Initialize.h"

int main()
{
	//驱动设置和中文路径设置
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* referImg = ".\\data\\reference.tif";
	const char* resultImg = ".\\data\\result_4677.tif";

	Initialize* refImg = new Initialize(referImg);
	Initialize* resImg = new Initialize(resultImg);


	EvaAccuracy(*refImg, *resImg);


	system("pause");
	return 0;
}


