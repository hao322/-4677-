#include "pch.h"
#include "DataIni.h"

int main()
{
    //驱动设置和中文路径设置
    GDALAllRegister();
    CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");

    const char* imgFile = ".\\data\\testImg.tif";

    // 打开图像获取数据
    DataIni* Img = new DataIni(imgFile);

    int Xsize = Img->getXsize();
    int Ysize = Img->getYsize();
    GDALDataset* poDataset = Img->getDataset();
    GDALDataType poDatatype = Img->getDatatype();

    unsigned char* imageData = new unsigned char[Xsize * Ysize * poDatatype];
    poDataset->GetRasterBand(1)->RasterIO(GF_Read, 0, 0x0, Xsize, Ysize, imageData, Xsize, Ysize, poDatatype, 0, 0);

    // KMeans分类和保存
    kmeansImage(imageData, *Img, 3, ".\\data\\result_4677.tif");

    // 清理
    delete[] imageData;
    GDALClose(poDataset);

    cout << "完成，请查看输出结果，保存于data目录下" << endl;

    system("pause");
    return -1;
}

