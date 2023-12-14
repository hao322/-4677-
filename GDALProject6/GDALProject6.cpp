#include "pch.h"
#include "Initialize.h"

int main()
{
	//驱动设置和中文路径设置
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* afterImg = ".\\data\\tsunami_after.tif";
	const char* beforeImg = ".\\data\\tsunami_before.tif";
	const char* resultImg = ".\\data\\result_4677_1.tif";

	Initialize* after = new Initialize(afterImg);
	Initialize* before = new Initialize(beforeImg);

	unsigned int* _Change = getChange(*after, *before);
	
	int Threshold = threshold(_Change, after->getXsize(), after->getYsize());

	unsigned short* result = segmentation(_Change, Threshold, after->getXsize(), after->getYsize());

	saveFile(result, resultImg, afterImg);

	system("pause");
	return 0;
}

