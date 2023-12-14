#include "pch.h"
#include <gdal_priv.h>
#include <iostream>
#include "Filtering.h"
using namespace std;


int main()
{
	//注册文件格式
	GDALAllRegister();
	//设置支持中文路径
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* imgpath = ".\\data\\nosieImg.tif";
	const char* savepath1 = ".\\data\\Noise_result_4677.tif";
	const char* savepath2 = ".\\data\\Noise_MedianResult.tif";
	Filtering* filter = new Filtering(imgpath);

	//保存——均值滤波
	//filter->SaveFile(savepath1, filter->Mean());
	//保存——中值滤波
	filter->SaveFile(savepath2, filter->Median());
	
	system("pause");
	return 0;

}