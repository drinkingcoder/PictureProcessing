#include <math.h>
#include<stdio.h>
#include<iostream>
#include<string>
#include<stdlib.h>
#include"../Defs.h"
#define PI 3.1415926

#ifndef _MATRIX_H_
#define _MATRIX_H_

using namespace std;

class Matrix
{
public:
	float** data;	
	int Row,Col;
	Matrix(int Row,int Col);
	Matrix(Matrix const& m);
	~Matrix();

	float* operator[](int k);
	friend Matrix operator+(Matrix& m1,Matrix& m2);
	friend Matrix operator-(Matrix& m1,Matrix& m2);
	friend Matrix operator*(Matrix& m1,Matrix& m2);

	void SwapRow(int i,int j);
	void TimesRow(int row,float times);
	void PlusRow(int i,float times,int j);
	void SwapCol(int i,int j);
	Matrix* SolveEquation();

	static float GaussianFunction(float sigma,float dis);

	static Matrix* UnitMatrix(int Row);
	static Matrix* Multiplication(Matrix* m1,Matrix* m2);
	static Matrix* Translate(Matrix* m,float x,float y);
	static Matrix* MirrorX(Matrix* m,float x);
	static Matrix* MirrorY(Matrix* m,float y);
	static Matrix* Scale(Matrix* m,float xcoeff,float ycoeff);
	static Matrix* Rotate(Matrix* m,float degree);
	static Matrix* Shear(Matrix* m,float x,float y);
	static Matrix* Laplacian(int Row);
	static Matrix* ones(int Row);
	static Matrix* GaussianFilter(float sigma);
	void printInfo(string s);
};



#endif
