#include"Matrix.h"
#include<stdio.h>
#include<iostream>

#define DEBUG

int main()
{
	Matrix* m = new Matrix(4,5);
	srand(time(NULL));
	for(int i=0;i<4;i++)
		for(int j=0;j<5;j++)
			m->data[i][j] = rand()%10;
	m->printInfo("\t");
	cout << "--------" << endl;
	m =	m->SolveEquation();
	m->printInfo("--");
}

