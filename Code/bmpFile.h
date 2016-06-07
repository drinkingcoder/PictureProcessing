#include<iostream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include "Matrix/Matrix.h"

#define BOXSIZE 100

using namespace std;

typedef struct bfHeader
{
	unsigned short bfType;			//2,4,2,2,4
	unsigned int bfSize;
	unsigned short bfReserved1;
	unsigned short bfReserved2;
	unsigned int bfOffBits;
} bfHeader;

typedef struct biHeader
{
	unsigned int biSize;			//4,4,4,2
	unsigned int biWidth;
	unsigned int biHeight;
	unsigned short biPlanes;
	unsigned short biBitCount;		//2,4,4,4
	unsigned int biCompression;
	unsigned int biSizeImage;
	unsigned int biXPixelsPerMeter;
	unsigned int biYPixelsPerMeter;	//4,4,4
	unsigned int biColorUsed;
	unsigned int biColorImportant;
} biHeader;	

typedef struct RGBE
{
	unsigned char Blue,Green,Red,Reserved;
} RGBE;

class bmpFile {
private:
	double transferRGBToGrey(int Blue,int Green,int Red);
	double calcBetweenVariance(int threshold,int maximal,int minimal,int startRow,int startCol,int calcWidth,int calcHeight);
	void assignBinaryValue(unsigned char** data,int deltaWidth,int deltaHeight,int threshold,int startRow,int startCol,int width,int height,int flag);
	void getMaxMin(int& maximal,int& minimal,int startRow,int startCol,int calcWidth,int calcHeight);
	int flagRow,flagCol;
	int sum[256];
	int amount[256];
	int row_amount[256];
	unsigned char** data;
	double max_between_variance;
	double wf,wb;
	double tot_amount,amountf,amountb;
    int minimal,maximal;
	double average_fore,average_back,sum_fore,sum_back;
	int chosen_threshold;
public:
	bfHeader bmpFileHeader;
	biHeader bmpInfoHeader;
	RGBE*	palette;
	unsigned char** bmpFileData;
	int** yuv;
	long nBytesPerLine;
	Matrix* transformMatrix;
	Matrix* reversTransformMatrix;

	bmpFile(char* fileName);
	~bmpFile();
	int exportToFile(char* fileName);
	float sqr(float);
	void printInfo();
	void transferToGreyBMP();
	void globalBinarization();
	void localBinarization(int width,int height,double constriction);
	void changeLuminance(int luminance);
	void dilation(unsigned char **m,int height,int width);
	void erotion(unsigned char **m,int height,int width);
	void open(unsigned char **m,int height,int width);
	void close(unsigned char **m,int height,int width);
	void complete();
	void extract(unsigned char **m,int height,int width);
	//the input height&width must be odd,and the checkpoint is at (height/2,width/2)
	//accept a matrix to operate
	void HistogramLinearEqualization();
	void HistogramLinearEqualizationWithLuminance();
	//res = start + k*logbase(sum[ori])
	//k = (finish - start )/logbase(tot)
	//restrict the data between 0~255
	void HistogramLogEqualization(float min_pixel,float max_pixel,float base);
	void HistogramLogEqualizationWithLuminance(float min_pixel,float max_pixel,float base);
	//no data restriction, if the data overflowed , it will be assign by 0 or 255;
	void HistogramLogEqualization2(float a,float b,float base);
	void HistogramLogEqualizationWithLuminance2(float a,float b,float base);
	void transferDataToYUV();
	void transferYUVToData();

	//Lab4 start
	void fill(int row,int col);
	void mapResultMatrixWithCanvasSize(float canvasScaleCoef);
	void mirrorX(int x);
	void mirrorY(int y);
	void translate(int x,int y);
	void scale(float x,float y);
	void rotate(float degree);
	void shear(float x,float y);

	//Lab5 start
	unsigned char clip(int k);
	unsigned char clip(float k);
	float Convolution(Matrix* m,int r,int c,int times,int k);
	//convolute the matrix with bmpFileData;point (r,c) is locate in the left up corner
	void LaplacianEnhancement();
	void MeanFilter(int length);
	void GaussianFilter(float sigma);

	//Lab5 start
	void generateBilateralFilter(Matrix* inputFilter,Matrix* resultFilter,float sigmar,int r,int c,int times,int k);
	void generateColorfulBilateralFilter(Matrix* inputFilter,Matrix* resultFilter,float sigmar,int r,int c);
	void BilateralFilter(float sigmas,float sigmar);
	void ColorfulBilateralFilter(float sigmas,float sigmar);
	float getVariance(Matrix* inputFilter,int r,int c,int times);
	void AutoBilateralFilter(float sigmas);
};

