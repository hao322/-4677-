// GDALProject2.cpp : 此文件包含 "main" 函数。程序执行将在此处开始并结束。
//

#include "pch.h"
#include "Data.h"
#include <iostream>
#include <gdal_priv.h>



int main()
{
	//注册文件格式
	GDALAllRegister();
	//设置支持中文路径
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* imgpath = ".\\data\\test.tif";
	const char* savepath = ".\\data\\vegetable_4677.tif";

	Data* data = new Data(imgpath);
	data->saveFile(savepath, data->getNDVI());

	system("pause");
	return 0;
}


