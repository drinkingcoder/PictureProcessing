#include "bmpFile.h"

#define PI 3.1415926
#define DEBUG

char fileName[100] = "/Users/drinkingcoder/Documents/university/Graph&Info/Lab6/autosigmar-20sigmas-3itr.bmp";
char Gaussian[100] = "/Users/drinkingcoder/Documents/university/Graph&Info/Lab5/Gaussian.bmp";
char Mean[100] = "/Users/drinkingcoder/Documents/university/Graph&Info/Lab5/Mean.bmp";
char fileout[100] = "/Users/drinkingcoder/Documents/university/Graph&Info/Lab6/gray.bmp";
char fileenhance[100] = "/Users/drinkingcoder/Documents/university/Graph&Info/Lab6/enhance.bmp";
int main()
{
	bmpFile* bf = new bmpFile(fileName);
//	bf->ColorfulBilateralFilter(10,10);
	bf->AutoBilateralFilter(20);
	bf->exportToFile(fileenhance);
	delete bf;

}
