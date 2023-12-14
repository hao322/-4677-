
#include "pch.h"
#include "HIS_ImgFusion.h"
#include <iostream>


int main()
{
	//驱动设置和中文路径设置
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//影像地址
	const char* Mul_Path = ".\\data\\多光谱.tif"; 
	const char* Pan_Path = ".\\data\\全色.tif";
    
	Tttt(Mul_Path, Pan_Path);


    system("pause");
    return 0;
}
