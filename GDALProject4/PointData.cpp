#include "pch.h"
#include "PointData.h"

PointData::PointData(const char* DataPath)
{
	this->Dataset = (GDALDataset*)GDALOpen(DataPath, GA_ReadOnly);
	this->Xsize = Dataset->GetRasterXSize();
	this->Ysize = Dataset->GetRasterYSize();
	this->Bandnum = Dataset->GetRasterCount();
	this->DataType = Dataset->GetRasterBand(1)->GetRasterDataType();
}
int PointData::getXsize()
{
	return this->Xsize;
}
int PointData::getYsize()
{
	return this->Ysize;
}
int PointData::getBandnum()
{
	return this->Bandnum;
}
PointData::~PointData()
{
	//cout << "程序执行完毕，请查看！" << endl;
}


//计算一次多项式
tuple<Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d, Eigen::Vector3d> PolynomialCal(const char* ptData)
{
	Eigen::VectorXd A1(18);
	Eigen::VectorXd A2(18);
	Eigen::Matrix<double, 18, 3> B;
	Eigen::VectorXd C1(18);
	Eigen::VectorXd C2(18);
	Eigen::Matrix<double, 18, 3> D;
	Eigen::Vector3d X1;
	Eigen::Vector3d X2;
	Eigen::Vector3d Y1;
	Eigen::Vector3d Y2;
	ifstream ifs;
	ifs.open(ptData, ios::in);
	if (ifs.is_open())
	{
		string buf;
		int i = 0;
		while (getline(ifs, buf))
		{
			double a = 0, b = 0, c = 0, d = 0;
			//类型转换，string To double
			istringstream str2double(buf); //自动以空格分隔
			str2double >> a >> b >> c >> d;
			A1[i] = a; A2[i] = b; B(i, 0) = 1; B(i, 1) = c; B(i, 2) = d;
			C1[i] = c; C2[i] = d; D(i, 0) = 1; D(i, 1) = a; D(i, 2) = b;
			i += 1;
		}
		//正运算系数
		X1 = (B.transpose() * B).inverse() * B.transpose() * A1;
		Y1 = (B.transpose() * B).inverse() * B.transpose() * A2;
		//反运算系数
		X2 = (D.transpose() * D).inverse() * D.transpose() * C1;
		Y2 = (D.transpose() * D).inverse() * D.transpose() * C2;

	}
	else
	{
		cout << "文件打开失败" << endl;
	}
	//cout << X1[0] << "  " << X1[1] << "  " << X1[2] << endl;
	//cout << Y1[0] << "  " << Y1[1] << "  " << Y1[2] << endl;

	//cout << X2[0] << "  " << X2[1] << "  " << X2[2] << endl;
	//cout << Y2[0] << "  " << Y2[1] << "  " << Y2[2] << endl;
	//关闭文件
	ifs.close();
	return make_tuple(X1, Y1, X2, Y2);
}

int ToMax(int a, int b, int c, int d)
{
	int max = 0;
	if (b > a && b > c && b > d)
		max = b;
	else if (c > a && c > b && c >d)
		max = c;
	else if (d > a && d > b && d > c)
		max = d;
	else
		max = a;
	return max;
}

//几何校正
void Geometric_Correction(PointData& Srcset, PointData& Corset, const char* pdata, const char* resultPath)
{
	//X1
	double a0 = get<0>(PolynomialCal(pdata))[0];
	double a1 = get<0>(PolynomialCal(pdata))[1];
	double a2 = get<0>(PolynomialCal(pdata))[2];
	double b0 = get<1>(PolynomialCal(pdata))[0];
	double b1 = get<1>(PolynomialCal(pdata))[1];
	double b2 = get<1>(PolynomialCal(pdata))[2];
	double a0_ = get<2>(PolynomialCal(pdata))[0];
	double a1_ = get<2>(PolynomialCal(pdata))[1];
	double a2_ = get<2>(PolynomialCal(pdata))[2];
	double b0_ = get<3>(PolynomialCal(pdata))[0];
	double b1_ = get<3>(PolynomialCal(pdata))[1];
	double b2_ = get<3>(PolynomialCal(pdata))[2];
	
	int msX = Srcset.Xsize;
	int msY = Srcset.Ysize;
	int mcX = Corset.Xsize;
	int mcY = Corset.Ysize;


	//利用校正多项式2计算待配准影像的四个角点的坐标
	int LTx = a0_ + a1_ * 0 + a2_ * 0;        int LTy = b0_ + b1_ * 0 + b2_ * 0;
	int RTx = a0_ + a1_ * mcX + a2_ * 0;      int RTy = b0_ + b1_ * mcX + b2_ * 0;
	int LBx = a0_+ a1_ * 0 + a2_ * mcY;      int LBy = b0_ + b1_ * 0 + b2_ * mcY;
	int RBx = a0_ + a1_ * mcX + a2_ * mcY;	   int RBy = b0_ + b1_ * mcX + b2_ * mcY;


	//求出构建X方向最大值、Y方向的最大值
	int maxX = 0; int maxY = 0;
	maxX = ToMax(LTx, RTx, LBx, RBx);
	maxY = ToMax(LTy, RTy, LBy, RBy);

	GDALDataType mDataType = Corset.Dataset->GetRasterBand(1)->GetRasterDataType();
	GDALRasterBand* mSrcBand = Srcset.Dataset->GetRasterBand(1);
	GDALRasterBand* mCorBand = Corset.Dataset->GetRasterBand(1);
	unsigned char* mSrcBand_ = (unsigned char*)CPLMalloc(msX * msY * mDataType);
	unsigned char* mCorBand_ = (unsigned char*)CPLMalloc(mcX * mcY * mDataType);
	//计算在参考影像上的位置，赋予灰度值
	unsigned char* resultbuf = (unsigned char*)CPLMalloc(sizeof(mDataType) * maxX * maxY);
	mSrcBand->RasterIO(GF_Read, 0, 0, msX, msY, mSrcBand_, msX, msY, mDataType, 0, 0);
	mCorBand->RasterIO(GF_Read, 0, 0, mcX, mcY, mCorBand_, mcX, mcY, mDataType, 0, 0);

	for (int i = 0; i < maxY; i++)
	{
		for (int j = 0; j < maxX; j++)
		{
			float x = a0 + a1 * j + a2 * i;
			float y = b0 + b1 * j + b2 * i;
			//不在图像的范围内，赋值为0
			if (x < 0 || x > msX || y < 0 || y > msY)
			{ 
				resultbuf[i * maxX + j] = 0;
			}
			else
			{
				int p1x = floor(x), p1y = floor(y);
				int p2x = ceil(x), p2y = floor(y);
				int p3x = floor(x), p3y = floor(y);
				int p4x = ceil(x), p4y = ceil(y);
				float d_t = y - floor(y);
				float d_b = ceil(y) - y;
				float d_l = x - floor(x);
				float d_r = ceil(x) - x;

				resultbuf[i * maxX + j] = mSrcBand_[p1y * msX + p1x] * d_b * d_r + mSrcBand_[p2y * msX + p2x] * d_b * d_l +
					mSrcBand_[p3y * msX + p3x] * d_t * d_r + mSrcBand_[p4y * msX + p4x] * d_t * d_l;
			}
		}
	}

	//获取驱动
	GDALDriver* Driver = GetGDALDriverManager()->GetDriverByName("GTiff");
	//定义波段排列顺序
	int Bandmap[1] = { 1 };

	//创建保存影像的数据集
	char** papszOptions = nullptr;
	papszOptions = CSLSetNameValue(papszOptions, "INTERLEAVE", "BAND");
	GDALDataset* imgDstDS = Driver->Create(resultPath , maxX, maxY, 1, mDataType, papszOptions);
	//若创建失败
	if (!imgDstDS)
	{
		cout << "Create file failed!" << endl;
	}
	else
	{
		//写入图像
		imgDstDS->GetRasterBand(1)->RasterIO(GF_Write, 0, 0, maxX, maxY, resultbuf, maxX, maxY, mDataType, 0, 0);
		cout << "Success，请查看文件！" << endl;
	}

	CPLFree(mCorBand_);
	CPLFree(mSrcBand_);
}

//NCC_match
void Match_Ncc(PointData& SrcData, PointData& RefData)
{
	const int r = 80;
	int msX = SrcData.Dataset->GetRasterXSize();
	int msY = SrcData.Dataset->GetRasterYSize();
	int mrX = RefData.Xsize;
	int mrY = RefData.Ysize;
	GDALDataType mDataType = SrcData.Dataset->GetRasterBand(1)->GetRasterDataType();
	GDALRasterBand* mSrcBand = SrcData.Dataset->GetRasterBand(1);
	GDALRasterBand* mRefBand = RefData.Dataset->GetRasterBand(1);
	unsigned char* _mSrcBand = (unsigned char*)CPLMalloc(msX * msY * mDataType);
	unsigned char* _mRefBand = (unsigned char*)CPLMalloc(mrX * mrY * mDataType);
	mSrcBand->RasterIO(GF_Read, 0, 0, msX, msY, _mSrcBand, msX, msY, mDataType, 0, 0);
	mRefBand->RasterIO(GF_Read, 0, 0, mrX, mrY, _mRefBand, mrX, mrY, mDataType, 0, 0);
	unsigned char* Temp = (unsigned char*)CPLMalloc((2 * r + 1) * (2 * r + 1) * mDataType);

	int sX = msX / 4, sY = msY / 4;
	int lx = sX - r, rx = sX + r;
	int ty = sY - r, by = sY + r;
	if (lx < 0 || rx > msX || ty < 0 || by > msY)
	{
		cout << "窗口过大！请退出重新设置窗口大小" << endl;
		system("pause");
	}
	else
	{
		for (int i = ty; i <= by; i++)
		{
			for (int j = lx; j <= rx; j++)
			{
				Temp[(i - ty) * (2 * r + 1) + (j - lx)] = _mSrcBand[i * msX + j];
			}
		}
	}

	//NCC匹配算法
	int rX, rY;
	int max_X = 0;
	int max_Y = 0;
	float max_ncc = -1;
	int mTmpX = 2 * r + 1;
	int mTmpY = 2 * r + 1;
	int k = 0;
	for (int i = 0; i < mrY - mTmpY + 1; i++)
		for (int j = 0; j < mrX - mTmpX + 1; j++)
		{
			float sum1 = 0;
			float sum2 = 0;
			float ave1 = 0;
			float ave2 = 0;
			for (int m = 0; m < mTmpY; m++)
				for (int n = 0; n < mTmpX; n++)
				{
					sum1 = sum1 + (float)_mRefBand[(i + m) * mrX + (j + n)];
					sum2 = sum2 + (float)Temp[m * mTmpX + n];
				}
			ave1 = sum1 / (mTmpX * mTmpY);
			ave2 = sum2 / (mTmpX * mTmpY);
			float sum_a_b = 0;
			float sum_a_2 = 0;
			float sum_b_2 = 0;
			for (int m = 0; m < mTmpY; m++)
				for (int n = 0; n < mTmpX; n++)
				{
					float a = (float)_mRefBand[(i + m) * mrX + (j + n)] - ave1;
					float b = (float)Temp[m * mTmpX + n] - ave2;
					sum_a_b += a * b;
					sum_a_2 += pow(a, 2);
					sum_b_2 += pow(b, 2);
				}
			float ncc = sum_a_b / sqrt(sum_a_2 * sum_b_2);
			if (ncc > max_ncc)
			{
				max_ncc = ncc;
				max_X = j;
				max_Y = i;
			}
			while (max_ncc > 0.5)
			{
				k += 1;
			}
		}
	cout << "找到同名点" << k << "个" << endl;
	rX = max_X + r;
	rY = max_Y + r;
	CPLFree(_mRefBand);
	CPLFree(_mSrcBand);
	CPLFree(Temp);

}

