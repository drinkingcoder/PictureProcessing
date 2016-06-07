#include "Matrix.h"

float sqr(float k)
{
	return k*k;
}

Matrix::Matrix(int Row,int Col)
{
	this->Row = Row;
	this->Col = Col;
	data = (float**)malloc(sizeof(float*)*Row);
	for(int i=0;i<Row;i++)
		data[i] = new float[Col];
	for(int i=0;i<Row;i++)
		for(int j=0;j<Col;j++)
			data[i][j] = 0;
}

Matrix* Matrix::UnitMatrix(int Row)
{
	Matrix* res = new Matrix(Row,Row);
	for(int i=0;i<Row;i++)
		res->data[i][i] = 1;
	return res;
}

Matrix::Matrix(Matrix const& m)
{
	Col = m.Col;
	Row = m.Row;
	data = (float**)malloc(sizeof(float*)*Row);
	for(int i=0;i<Row;i++)
		data[i] = new float[Col];
	for(int i=0;i<Row;i++)
		for(int j=0;j<Col;j++)
			data[i][j] = m.data[i][j];
}

Matrix::~Matrix()
{
	for(int i=0;i<Row;i++)
		free(data[i]);
	free(data);
}

float* Matrix::operator[](int k)
{
	return data[k];
}

Matrix* Matrix::Multiplication(Matrix* m1,Matrix* m2)
{
	if(m1->Col != m2->Row)
	{
		#ifdef DEBUG
			cout << "Matrix::Multiplication" << endl;
			cout << "	m1->Col != m2->Row" << endl;
			cout << "	m1.Col = " << m1->Col << endl;
			cout << "	m2.Row = " << m2->Row << endl;
		#endif
	  	return NULL;
	}
	Matrix* res = new Matrix(m1->Row,m2->Col);
	for(int i=0;i<m1->Row;i++)
		for(int k=0;k<m1->Col;k++)
			for(int j=0;j<m2->Col;j++)
				res->data[i][j] += m1->data[i][k]*m2->data[k][j];
	delete m1;
	delete m2;
	return res;
}

Matrix* Matrix::Translate(Matrix* m,float x,float y)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::Translate" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* translateMatrix = new Matrix(3,3);
	translateMatrix->data[0][0] = 1;
	translateMatrix->data[1][1] = 1;
	translateMatrix->data[2][2] = 1;
	translateMatrix->data[0][2] = x;
	translateMatrix->data[1][2] = y;
	Matrix* resultMatrix = Matrix::Multiplication(translateMatrix,m);
	return resultMatrix;
}

Matrix* Matrix::MirrorX(Matrix* m,float x)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::MirrorX" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* resultMatrix;
	Matrix* mirrorMatrix = new Matrix(3,3);
	mirrorMatrix->data[0][0] = -1;
	mirrorMatrix->data[1][1] = 1;
	mirrorMatrix->data[2][2] = 1;

	mirrorMatrix = Matrix::Translate(mirrorMatrix,2*x,0);
	resultMatrix = Matrix::Multiplication(mirrorMatrix,m);
	return resultMatrix;
}

Matrix* Matrix::MirrorY(Matrix* m,float y)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::MirrorY" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* mirrorMatrix = new Matrix(3,3);
	Matrix* resultMatrix;

	mirrorMatrix->data[0][0] = 1;
	mirrorMatrix->data[1][1] = -1;
	mirrorMatrix->data[2][2] = 1;

	mirrorMatrix = Matrix::Translate(mirrorMatrix,0,2*y);
	resultMatrix = Matrix::Multiplication(mirrorMatrix,m);

	return resultMatrix;
}

Matrix* Matrix::Scale(Matrix* m,float xcoeff,float ycoeff)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::Scale" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* scaleMatrix = new Matrix(3,3);
	Matrix* resultMatrix;

	scaleMatrix->data[0][0] = xcoeff;
	scaleMatrix->data[1][1] = ycoeff;
	scaleMatrix->data[2][2] = 1;

	resultMatrix = Matrix::Multiplication(scaleMatrix,m);

	return resultMatrix;
}

Matrix* Matrix::Rotate(Matrix* m,float degree)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::Rotate" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* rotateMatrix = new Matrix(3,3);
	Matrix* resultMatrix;

	rotateMatrix->data[0][0] = rotateMatrix->data[1][1] = cos(degree);
	rotateMatrix->data[0][1] = -sin(degree);
	rotateMatrix->data[1][0] = sin(degree);
	rotateMatrix->data[2][2] = 1;

	resultMatrix = Matrix::Multiplication(rotateMatrix,m);

	return resultMatrix;
}

Matrix* Matrix::Shear(Matrix* m,float x,float y)
{
	if(m->Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::Shear" << endl;
			cout << "	m.Row = " << m->Row << endl;
		#endif
		return NULL;
	}
	Matrix* rotateMatrix = new Matrix(3,3);
	Matrix* resultMatrix;

	rotateMatrix->data[0][0] = rotateMatrix->data[1][1] = 1;
	rotateMatrix->data[0][1] = x;
	rotateMatrix->data[1][0] = y;
	rotateMatrix->data[2][2] = 1;

	resultMatrix = Matrix::Multiplication(rotateMatrix,m);

	return resultMatrix;
}

Matrix* Matrix::Laplacian(int Row)
{
	if(Row !=3 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::Laplacian" << endl;
			cout << "	Row = " << Row << endl;
		#endif
		return NULL;
	}
	Matrix* resultMatrix = new Matrix(Row,Row);
	resultMatrix->data[0][0] = resultMatrix->data[0][2] = -1;
	resultMatrix->data[2][0] = resultMatrix->data[2][2] = -1;
	resultMatrix->data[0][1] = resultMatrix->data[2][1] = -1;
	resultMatrix->data[1][0] = resultMatrix->data[1][2] = -1;
	resultMatrix->data[1][1] = 8;
	return resultMatrix;
}

Matrix* Matrix::ones(int Row)
{
	if(Row <= 0 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::ones" << endl;
			cout << "	Row = " << Row << endl;
		#endif
		return NULL;
	}
	Matrix* resultMatrix = new Matrix(Row,Row);
	for(int i=0;i<Row;i++)
		for(int j=0;j<Row;j++)
			resultMatrix->data[i][j] = 1;	
	return resultMatrix;
}

float Matrix::GaussianFunction(float sigma,float dis)
{
	return	1/(sigma*sqrt(2*PI))*exp(-sqr(dis)/(2*sqr(sigma)));
}

Matrix* Matrix::GaussianFilter(float sigma)
{
	if(sigma <= 0 ) 
	{
		#ifdef DEBUG
			cout << "Matrix::GaussianFilter" << endl;
			cout << "	Sigma = " << sigma << endl;
		#endif
		return NULL;
	}
	Matrix* resultMatrix = new Matrix(int(sigma+0.5),int(sigma+0.5));
	float center = resultMatrix->Row*1.0/2;
	float dis;
	for(int i=0;i<resultMatrix->Row;i++)
		for(int j=0;j<resultMatrix->Col;j++)
		{
			dis = sqrt(sqr(center-i-0.5)+sqr(center-j-0.5));
			resultMatrix->data[i][j] =  GaussianFunction(sigma,dis);
		}
	return resultMatrix;
}

Matrix* Matrix::SolveEquation()
{
	if(Row != Col-1)
	{
		#ifdef DEBUG
			cout<< "Matrix::SolveEquation" << endl;
			cout << "	m.Row = " << Row << endl;
			cout << "	m.Col = " << Col << endl;
		#endif
		return NULL;
	}
	int i,j,k;
	float times;
	Matrix* m = new Matrix(Row,1);
	for(i=0;i<Row;i++)
	{
		for(j=i;j<Row;j++)
			if(data[j][i]!=0)
			{
				SwapRow(i,j);
				break;
			}
		if(data[i][i] == 0) continue;
		for(int j=i+1;j<Row;j++)
			if(data[j][i]!=0)
				PlusRow(j,-data[j][i]/data[i][i],i);
	}
	for(int i=Row-1;i>=0;i--)
	{
		if(data[i][i] == 0) 
		{
			delete m;
			return NULL;
		}
		data[i][Col-1]/=data[i][i];
		data[i][i] = 1;
		for(int j=0;j<i;j++)
			PlusRow(j,-data[j][i],i);
		m->data[i][0] = data[i][Col-1];
	}
#ifdef DEBUG
	cout << "Solved Ans = " << endl;
	m->printInfo("\t");
#endif 
	return m;
}

void Matrix::SwapRow(int i,int j)
{
	if(i<0 || i>=Row)
	{
		#ifdef DEBUG
			cout<< "Matrix::SwapRow" << endl;
			cout << "	i = " << i << endl;
		#endif
		return;
	}
	if(j<0 || j>=Row)
	{
		#ifdef DEBUG
			cout<< "Matrix::SwapRow" << endl;
			cout << "	j = " << j << endl;
		#endif
		return;
	}
	float tmp;
	for(int k=0;k<Col;k++)
	{
		tmp = data[i][k];
		data[i][k] = data[j][k];
		data[j][k] = tmp;
	}
}

void Matrix::TimesRow(int row,float times)
{
	if(row<0 || row>=Row)
	{
		#ifdef DEBUG
			cout<< "Matrix::TimesRow" << endl;
			cout << "	row = " << row << endl;
		#endif
		return;
	}
	for(int i=0;i<Col;i++)
		data[row][i]*=times;	
	return;
}

void Matrix::PlusRow(int i,float times,int j)
{
	if(i<0 || i>=Row)
	{
		#ifdef DEBUG
			cout<< "Matrix::PlusRow" << endl;
			cout << "	i = " << i << endl;
		#endif
		return;
	}
	if(j<0 || j>=Row)
	{
		#ifdef DEBUG
			cout<< "Matrix::PlusRow" << endl;
			cout << "	j = " << j << endl;
		#endif
		return;
	}
	for(int k=0;k<Col;k++)
		data[i][k]+=times*data[j][k];
}

void Matrix::SwapCol(int i,int j)
{
	if(i<0 || i>=Col)
	{
		#ifdef DEBUG
			cout<< "Matrix::SwapCol" << endl;
			cout << "	i = " << i << endl;
		#endif
		return;
	}
	if(j<0 || j>=Col)
	{
		#ifdef DEBUG
			cout<< "Matrix::SwapCol" << endl;
			cout << "	j = " << j << endl;
		#endif
		return;
	}
	float tmp;
	for(int k=0;k<Row;k++)
	{
		tmp = data[k][i];
		data[k][i] = data[k][j];
		data[k][j] = tmp;
	}
}

Matrix operator*(Matrix& m1,Matrix& m2)
{
	Matrix res(m1.Row,m2.Col);
	if(m1.Col != m2.Row)
	{
		#ifdef DEBUG
			cout << "Matrix::operator* " << endl;
			cout << "	m1.Col != m2.Row" << endl;
			cout << "	m1.Col = " << m1.Col << endl;
			cout << "	m2.Row = " << m2.Row << endl;
		#endif
		return res;
	}
	for(int i=0;i<m1.Row;i++)
		for(int k=0;k<m1.Col;k++)
			for(int j=0;j<m2.Col;j++)
				res[i][j] += m1[i][k]*m2[k][j];
	return res;
}

Matrix operator-(Matrix& m1,Matrix& m2)
{
	Matrix res(m1.Row,m1.Col);
	if(m1.Row != m2.Row || m1.Col != m2.Col)
	{
		#ifdef DEBUG
			cout << "Matrix::operator- " << endl;
			cout << "	m1.Row = " << m1.Row << endl;
			cout << "	m1.Col = " << m1.Col << endl;
			cout << "	m2.Row = " << m2.Row << endl;
			cout << "	m2.Col = " << m2.Col << endl;
		#endif
		return res;
	}
	for(int i=0;i<m1.Row;i++)
		for(int j=0;j<m1.Col;j++)
			res[i][j] = m1[i][j]-m2[i][j];
	return res;
}

Matrix operator+(Matrix& m1,Matrix& m2)
{
	Matrix res(m1.Row,m1.Col);
	if(m1.Row != m2.Row || m1.Col != m2.Col)
	{
		#ifdef DEBUG
			cout << "Matrix::operator+ " << endl;
			cout << "	m1.Row = " << m1.Row << endl;
			cout << "	m1.Col = " << m1.Col << endl;
			cout << "	m2.Row = " << m2.Row << endl;
			cout << "	m2.Col = " << m2.Col << endl;
		#endif
		return res;
	}
	for(int i=0;i<m1.Row;i++)
		for(int j=0;j<m1.Col;j++)
			res[i][j] = m1[i][j]+m2[i][j];
	return res;
}

void Matrix::printInfo(string s)
{
	for(int i=0;i<Row;i++)
	{
		cout << s;
		for(int j=0;j<Col;j++)
			cout << data[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}
