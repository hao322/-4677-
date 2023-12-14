
/*
根据模板影像，在一幅图像中寻找其最匹配位置，并输出其位置（行列号）
模板影像：".\\data\\template.tif"
图像:".\\data\\img.tif"
要求：使用基于NCC的图像匹配函数命名为 match_ncc
*/

#include "pch.h"
#include <iostream>
#include "Data.h"
#include <gdal_priv.h>

using namespace std;

int main()
{
	//注册文件格式
	GDALAllRegister();
	//设置支持中文路径
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//设置栅格文件路径
	const char* imgpath = ".\\data\\img.tif";
	const char* templatepath = ".\\data\\template.tif";

	Data* ImgData = new Data(imgpath);
	Data* TemData = new Data(templatepath);
	int Row, Col;
	Row = ImgData->match_ncc(*TemData).first;
	Col = ImgData->match_ncc(*TemData).second;
	cout << "匹配成功！" << endl;
	cout << "最佳匹配位置：" << Row << "," << Col << endl;

    system("pause");
    return 0;
}


