#include "bmpFile.h"
#include "Defs.h"

bmpFile::bmpFile(char* fileName)
{
#ifdef DEBUG
	cout << "bmpFile \'" << fileName << "\' creation start.." << endl;
#endif
	FILE* fp = fopen(fileName,"r");

	if( fp==NULL ) printf("Can't Open It!\n");
//bmpfileheader
	fread(&bmpFileHeader.bfType,2,1,fp);
	fread(&bmpFileHeader.bfSize,4,1,fp);
	fread(&bmpFileHeader.bfReserved1,2,1,fp);
	fread(&bmpFileHeader.bfReserved2,2,1,fp);
	fread(&bmpFileHeader.bfOffBits,4,1,fp);

//bmpinfoheader
	fread(&bmpInfoHeader.biSize,4,1,fp);
	fread(&bmpInfoHeader.biWidth,4,1,fp);
	fread(&bmpInfoHeader.biHeight,4,1,fp);
	fread(&bmpInfoHeader.biPlanes,2,1,fp);

	fread(&bmpInfoHeader.biBitCount,2,1,fp);
	fread(&bmpInfoHeader.biCompression,4,1,fp);
	fread(&bmpInfoHeader.biSizeImage,4,1,fp);
	fread(&bmpInfoHeader.biXPixelsPerMeter,4,1,fp);

	fread(&bmpInfoHeader.biYPixelsPerMeter,4,1,fp);
	fread(&bmpInfoHeader.biColorUsed,4,1,fp);
	fread(&bmpInfoHeader.biColorImportant,4,1,fp);

//palette
	if(bmpInfoHeader.biBitCount == 24) palette = NULL;
	else {
		palette =(RGBE*)malloc(sizeof(RGBE)*bmpInfoHeader.biBitCount);
		for( int i=0;i<bmpInfoHeader.biBitCount;i++)
			fread(&palette[i],4,1,fp);
	}
//nbytesperline
	switch(bmpInfoHeader.biBitCount)
	{
		case 1:nBytesPerLine=(bmpInfoHeader.biWidth/8+3)/4*4;break;
		case 4:nBytesPerLine=(bmpInfoHeader.biWidth/2+3)/4*4;break;
		case 8:nBytesPerLine=(bmpInfoHeader.biWidth+3)/4*4;break;
		case 16:nBytesPerLine=(bmpInfoHeader.biWidth*2+3)/4*4;break;
		case 24:nBytesPerLine=(bmpInfoHeader.biWidth*3+3)/4*4;break;
		default:printf("biBitCount is not proper!\n");break;
	}
//bmpFileData
	bmpFileData = (unsigned char**)malloc(sizeof(unsigned char*)*bmpInfoHeader.biHeight);

	data = (unsigned char**)malloc(sizeof(unsigned char*)*bmpInfoHeader.biHeight);
	yuv = (int**)malloc(sizeof(int*)*bmpInfoHeader.biHeight);
	transformMatrix = new Matrix(3,3);
	reversTransformMatrix = new Matrix(3,3);
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
	{
		bmpFileData[i]=(unsigned char*)malloc(nBytesPerLine);
		data[i]=(unsigned char*)malloc(nBytesPerLine);
		yuv[i] = (int*)malloc(nBytesPerLine*sizeof(int));
		for(int j=0;j<nBytesPerLine;j++)
			fread(&bmpFileData[i][j],1,1,fp);
	}
	fclose(fp);
	transformMatrix->data[0][0] = transformMatrix->data[1][1] = transformMatrix->data[2][2] = 1;
	reversTransformMatrix->data[0][0] = reversTransformMatrix->data[1][1] = reversTransformMatrix->data[2][2] = 1;
#ifdef DEBUG
	cout << "bmpFile \'" << fileName << "\' creation finished.." << endl;
#endif
}

double bmpFile::transferRGBToGrey(int Blue,int Red,int Green)
{
	return double(Blue*0.114+Green*0.587+Red*0.299);
}

void bmpFile::transferToGreyBMP()
{
	if(bmpInfoHeader.biBitCount==24)
	{
		bmpInfoHeader.biBitCount = 8;
		nBytesPerLine=(bmpInfoHeader.biWidth+3)/4*4;
		bmpInfoHeader.biSizeImage = nBytesPerLine*bmpInfoHeader.biHeight;
		bmpFileHeader.bfOffBits = 54+1024;
		bmpFileHeader.bfSize = bmpFileHeader.bfOffBits + bmpInfoHeader.biSizeImage;
		palette = (RGBE*)malloc(1024);
		for(int i=0;i<256;i++)
		{
			palette[i].Blue = palette[i].Red = palette[i].Green = i;
			palette[i].Reserved = 0;
		}
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j] = (int)transferRGBToGrey(bmpFileData[i][j*3],bmpFileData[i][j*3+1],bmpFileData[i][j*3+2]);
		return;
	} else 
		for(int i=0;i<(1<<bmpInfoHeader.biBitCount);i++)
			palette[i].Blue = palette[i].Green = palette[i].Red = transferRGBToGrey(palette[i].Blue,palette[i].Green,palette[i].Red);
}

void bmpFile::transferDataToYUV()
{
	int y,u,v;
	int r,g,b;
	if(bmpInfoHeader.biBitCount==24)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				b = bmpFileData[i][j*3];
				g = bmpFileData[i][j*3+1];
				r = bmpFileData[i][j*3+2];
				y = ( (  66 * r + 129 * g +  25 * b + 128) >> 8) ;
				u = ( ( -38 * r -  74 * g + 112 * b + 128) >> 8) ;
				v = ( ( 112 * r -  94 * g -  18 * b + 128) >> 8) ;
				yuv[i][j*3] = y;
				yuv[i][j*3+1] = u;
				yuv[i][j*3+2] = v;
			}
		return;
	}
}

void bmpFile::transferYUVToData()
{
	int y,u,v;
	int r,g,b;
	if(bmpInfoHeader.biBitCount==24)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				y = yuv[i][j*3] ;
				u = yuv[i][j*3+1];
				v = yuv[i][j*3+2];
				r = (int)(1.164383*y + 1.596027*v);
				g = (int)(1.164383*y - 0.391762*u - 0.812968*v);
				b = (int)(1.164383*y + 2.017232*u);
				if(b>255) b = 255;
				if(g>255) g = 255;
				if(r>255) r = 255;
				if(r<0) r = 0;
				if(g<0) g = 0;
				if(b<0) b = 0;
				bmpFileData[i][j*3] = (unsigned char)b;
				bmpFileData[i][j*3+1] = (unsigned char)g;
				bmpFileData[i][j*3+2] = (unsigned char)r;
			}
		return;
	}
}

void bmpFile::changeLuminance(int luminance)
{
#ifdef DEBUG
	cout << "changeLuminance start.." << endl;
#endif
	if(bmpInfoHeader.biBitCount==24)
	{
		transferDataToYUV();
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				yuv[i][j*3]+=luminance;
		transferYUVToData();
#ifdef DEBUG
		cout << "changeLuminance 24bits finished.." << endl;
#endif
		return;
	}
#ifdef DEBUG
	cout << "bmpFile not 24bits"<<endl;
#endif
}

bmpFile::~bmpFile()
{
#ifdef DEBUG
	cout << "destruct bmp start.." << endl;
	cout << "width = " << bmpInfoHeader.biWidth << endl;
	cout << "height = " << bmpInfoHeader.biHeight << endl;
#endif
	/*
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		free(data[i]);
	free(data);
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		free(bmpFileData[i]);
	free(bmpFileData);
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		free(yuv[i]);
	free(yuv);
	*/
	if(palette!=NULL)free(palette);
#ifdef DEBUG
	cout << "destruct bmp finished.." << endl;
#endif
}

int bmpFile::exportToFile(char* fileName)
{
#ifdef DEBUG
	cout << "exportToFile \'" << fileName << "\' start.." << endl;
#endif
	FILE* fp=fopen(fileName,"w");

	fwrite(&bmpFileHeader.bfType,2,1,fp);
	fwrite(&bmpFileHeader.bfSize,4,1,fp);
	fwrite(&bmpFileHeader.bfReserved1,2,1,fp);
	fwrite(&bmpFileHeader.bfReserved2,2,1,fp);
	fwrite(&bmpFileHeader.bfOffBits,4,1,fp);

	fwrite(&bmpInfoHeader.biSize,4,1,fp);
	fwrite(&bmpInfoHeader.biWidth,4,1,fp);
	fwrite(&bmpInfoHeader.biHeight,4,1,fp);
	fwrite(&bmpInfoHeader.biPlanes,2,1,fp);

	fwrite(&bmpInfoHeader.biBitCount,2,1,fp);
	fwrite(&bmpInfoHeader.biCompression,4,1,fp);
	fwrite(&bmpInfoHeader.biSizeImage,4,1,fp);
	fwrite(&bmpInfoHeader.biXPixelsPerMeter,4,1,fp);

	fwrite(&bmpInfoHeader.biYPixelsPerMeter,4,1,fp);
	fwrite(&bmpInfoHeader.biColorUsed,4,1,fp);
	fwrite(&bmpInfoHeader.biColorImportant,4,1,fp);
	switch(bmpInfoHeader.biBitCount)
	{
		case 1:fwrite(palette,2*4,1,fp);break;
		case 4:fwrite(palette,16*4,1,fp);break;
		case 8:fwrite(palette,256*4,1,fp);break;
		default:break;
	}
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		for(int j=0;j<nBytesPerLine;j++)
			fwrite(&bmpFileData[i][j],1,1,fp);
	fclose(fp);
#ifdef DEBUG
	cout << "exportToFile \'" << fileName << "\' finished.." << endl;
#endif
	return 1;
}

void bmpFile::printInfo()
{
	printf("bmpFileHeader:\n");
	printf("\tbfType = \t%u\n",bmpFileHeader.bfType);
	printf("\tbfSize = \t%u\n",bmpFileHeader.bfSize);
	printf("\tbfReserved1 = \t%u\n",bmpFileHeader.bfReserved1);
	printf("\tbfReserved2 = \t%u\n",bmpFileHeader.bfReserved2);
	printf("\tbfOffBits = \t%u\n",bmpFileHeader.bfOffBits);

	printf("bmpFileInfomation:\n");

	printf("\tbiSize = \t%u\n",bmpInfoHeader.biSize);
	printf("\tbiWidth = \t%d\n",bmpInfoHeader.biWidth);
	printf("\tbmpInfoHeadereight = \t%d\n",bmpInfoHeader.biHeight);
	printf("\tbiPlanes = \t%d\n",bmpInfoHeader.biPlanes);
	
	printf("\tbibitCount = \t%u\n",bmpInfoHeader.biBitCount);
	printf("\tbiCompression = \t%u\n",bmpInfoHeader.biCompression);
	printf("\tbiSizeImage = \t%u\n",bmpInfoHeader.biSizeImage);
	printf("\tbiXPixelsPerMeter = \t%u\n",bmpInfoHeader.biXPixelsPerMeter);

	printf("\tbiYPixelsPerMeter = \t%u\n",bmpInfoHeader.biYPixelsPerMeter);
	printf("\tbiColorUser = \t%u\n",bmpInfoHeader.biColorUsed);
	printf("\tbiColorImportant = \t%u\n",bmpInfoHeader.biColorImportant);

	if(palette == NULL) printf("Palette is NULL!\n");
	printf("size of imagedata = %d\n",nBytesPerLine*bmpInfoHeader.biHeight);
}

float bmpFile::sqr(float k)
{
	return k*k;
}

double bmpFile::calcBetweenVariance(int threshold,int maximal,int minimal,int startRow,int startCol,int calcWidth,int calcHeight)
{
	int finishRow = startRow+calcHeight;
	int finishCol = startCol+calcWidth;
	if(finishRow>bmpInfoHeader.biHeight) finishRow = bmpInfoHeader.biHeight;
	if(finishCol>bmpInfoHeader.biWidth)  finishCol = bmpInfoHeader.biWidth;
	double sum_fore = sum[threshold],sum_back = sum[maximal]-sum[threshold];
	double amount_fore = amount[threshold],amount_back = amount[maximal]-amount[threshold];
	double average_fore = sum_fore/amount_fore;
	double average_back = sum_back/amount_back;
	double amount_tot = amount[maximal];
	double between_variance = (amount_fore*amount_back)/sqr(amount_tot)*sqr(average_fore-average_back);
	return between_variance;
}

void bmpFile::getMaxMin(int& maximal,int& minimal,int startRow,int startCol,int calcWidth,int calcHeight)
{
	int finishRow = startRow+calcHeight;
	int finishCol = startCol+calcWidth;
	if(finishRow>bmpInfoHeader.biHeight) finishRow = bmpInfoHeader.biHeight;
	if(finishCol>bmpInfoHeader.biWidth)  finishCol = bmpInfoHeader.biWidth;
	maximal = 0;
	minimal = 255;
	for(int i=0;i<256;i++)
		sum[i] = amount[i] = 0;
	for(int	row = startRow; row < finishRow; row++)
		for(int col = startCol; col < finishCol; col++)
		{
			if(minimal>bmpFileData[row][col]) minimal = bmpFileData[row][col];
			if(maximal<bmpFileData[row][col]) maximal = bmpFileData[row][col];
			amount[bmpFileData[row][col]]++;
		}
	for(int i=1;i<=maximal;i++)
	{
		sum[i] = amount[i]*i+sum[i-1];
		amount[i]+=amount[i-1];
	}
}

void bmpFile::assignBinaryValue(unsigned char** data,int deltaWidth,int deltaHeight,int threshold,int startRow,int startCol,int width,int height,int flag)
{
	int finishRow = startRow+deltaHeight;
	int finishCol = startCol+deltaWidth;
	if(flagRow) finishRow = startRow+height;
	if(flagCol) finishCol = startCol+width;
	for(int row = startRow; row<finishRow;row++)
		for(int col = startCol; col<finishCol;col++)
		{
			if(flag) {
				data[row][col] = 255;
				continue;
			}
			if(bmpFileData[row][col]<threshold)
				data[row][col]= 0;
			else 
				data[row][col]= 255;
		}
}

void bmpFile::globalBinarization()
{
#ifdef DEBUG
	cout << "globalBinarization start.." << endl;
#endif
	transferToGreyBMP();
	
	int maximal,minimal;
	getMaxMin(maximal,minimal,0,0,bmpInfoHeader.biWidth,bmpInfoHeader.biHeight);
	double max_between_variance = -1;
	int chosen_threshold = minimal;
	for(int threshold = minimal; threshold < maximal ; threshold++)
	{
		double between_variance = calcBetweenVariance(threshold,maximal,minimal,0,0,bmpInfoHeader.biWidth,bmpInfoHeader.biHeight);
		if(between_variance >= max_between_variance)
		{
			max_between_variance = between_variance;
			chosen_threshold = threshold;
		}
	}
	assignBinaryValue(bmpFileData,bmpInfoHeader.biWidth,bmpInfoHeader.biHeight,chosen_threshold,0,0,bmpInfoHeader.biWidth,bmpInfoHeader.biHeight,0);
#ifdef DEBUG
	cout << "globalBinarization finished.." << endl;
#endif
}

void bmpFile::localBinarization(int width,int height,double constriction)
{
#ifdef DEBUG
	cout << "localBinarization start.." << endl;
#endif
	transferToGreyBMP();
	for(int row = 0; row < bmpInfoHeader.biHeight; row++)
		memcpy(data[row],bmpFileData[row],nBytesPerLine);
	int delta = 1;
    
    memset(amount,0,sizeof(amount));
    for(int	row = 0; row < height; row++)
        for(int col = 0; col < width; col++)
            amount[bmpFileData[row][col]]++;
    memcpy(row_amount,amount,sizeof(amount));

	flagRow = 0;
	for(int	row = 0; row < bmpInfoHeader.biHeight; row+=delta)
	{
		flagCol = 0;
		if(row >= bmpInfoHeader.biHeight-height)
		{
			flagRow = 1;
			row = bmpInfoHeader.biHeight-height;
		}
		for(int col = 0; col < bmpInfoHeader.biWidth; col+=delta)
		{
			if(col>=bmpInfoHeader.biWidth-width)
			{
				flagCol = 1;
				col = bmpInfoHeader.biWidth-width;
			}
            tot_amount =height*width;
            sum_back = 0;
            minimal = 255;
            maximal = 0;
            for(int i=0; i<256; i++)
            {
                sum_back+= row_amount[i]*i;
                if(row_amount[i]>0)
                {
                    if(minimal>i) minimal = i;
                    if(maximal<i) maximal = i;
                }
            }
            amountf = row_amount[minimal];
            wf = amountf/tot_amount;
            average_fore = minimal;
            chosen_threshold = minimal;
            sum_fore = minimal*row_amount[minimal];
            
            amountb = tot_amount - amountf;
            wb = amountb/tot_amount;
            sum_back -= minimal*row_amount[minimal];
            average_back = sum_back/amountb;
            max_between_variance = sqr(average_fore-average_back)*wb*wf;
			int chosen_threshold = minimal;
            for(int threshold = minimal+1; threshold < maximal ; threshold++)
            {
                amountf+= row_amount[threshold];
                amountb-= row_amount[threshold];
                wf = amountf/tot_amount;
                wb = amountb/tot_amount;
                sum_fore += row_amount[threshold]*threshold;
                sum_back -= row_amount[threshold]*threshold;
                average_fore = sum_fore/amountf;
                average_back = sum_back/amountb;
                double between_variance = sqr(average_fore-average_back)*wb*wf;
                if(between_variance > max_between_variance)
                {
                    max_between_variance = between_variance;
                    chosen_threshold = threshold;
                }
            }
			int flag = 0;
			if(max_between_variance<constriction) flag =1;
			assignBinaryValue(data,delta,delta,chosen_threshold,row,col,width,height,flag);
			if(flagCol) break;
            for(int i=row;i<row+height;i++)
            {
                row_amount[bmpFileData[i][col]]--;
                if((col+width)<bmpInfoHeader.biWidth)
                    row_amount[bmpFileData[i][col+width]]++;
            }
		}
		if(flagRow) break;
        for(int j=0;j<width;j++)
        {
            amount[bmpFileData[row][j]]--;
            if((row+height)<bmpInfoHeader.biHeight) amount[bmpFileData[row+height][j]]++;
        }
        memcpy(row_amount,amount,sizeof(amount));
	}
	for(int row = 0; row < bmpInfoHeader.biHeight; row++)
		memcpy(bmpFileData[row],data[row],nBytesPerLine);
#ifdef DEBUG
	cout << "localBinarization finished.." << endl;
#endif
}

void bmpFile::erotion(unsigned char **m,int height,int width)
{
#ifdef DEBUG
	cout << "erotion start.." << endl;
#endif
	int half_width = width/2;
	int half_height = height/2;
	int center_row,center_col;
	for(int row = 0; row<bmpInfoHeader.biHeight; row++)
		memset(data[row],0,nBytesPerLine);
	int hit;
	for(int row = 0; row <= bmpInfoHeader.biHeight-height; row++)
		for(int col = 0; col <= bmpInfoHeader.biWidth-width; col++)
		{
			hit = 1;
			center_row = row+half_height;
			center_col = col+half_width;	
			for(int i=0;i<height;i++)
			{
				for(int j=0;j<width;j++)
					if(m[i][j]^bmpFileData[row+i][col+j])
					{
						hit = 0;
						break;
					}
				if(hit == 0) break;
			}
			if(hit) data[row][col] = 255;
		}
	for(int row = 0; row < bmpInfoHeader.biHeight; row++)
		memcpy(bmpFileData[row],data[row],nBytesPerLine);
#ifdef DEBUG
	cout << "erotion finished.." << endl;
#endif
}

void bmpFile::dilation(unsigned char **m,int height,int width)
{
#ifdef DEBUG
	cout << "dilation start.." << endl;
#endif
	int half_width = width/2;
	int half_height = height/2;
	int center_row,center_col;
	for(int row = 0; row<bmpInfoHeader.biHeight; row++)
		memset(data[row],0,nBytesPerLine);
	for(int row = 0; row <= bmpInfoHeader.biHeight-height; row++)
		for(int col = 0; col <= bmpInfoHeader.biWidth-width; col++)
		{
			center_row = row+half_height;
			center_col = col+half_width;	
			if(bmpFileData[center_row][center_col] == 255)
				for(int i=0;i<height;i++)
					for(int j=0;j<width;j++)
						data[row+i][col+j] = 255;
		}
	for(int row = 0; row < bmpInfoHeader.biHeight; row++)
		memcpy(bmpFileData[row],data[row],nBytesPerLine);
#ifdef DEBUG
	cout << "dilation finished.." << endl;
#endif
}

void bmpFile::complete()
{
	for(int row=0; row<bmpInfoHeader.biHeight; row++)
		for(int col=0; col<bmpInfoHeader.biWidth; col++)
			if(bmpFileData[row][col] == 255) data[row][col] =0;
			else data[row][col] = 255;
	for(int row = 0; row < bmpInfoHeader.biHeight; row++)
		memcpy(bmpFileData[row],data[row],nBytesPerLine);
}

void bmpFile::close(unsigned char **m,int height,int width)
{
	erotion(m,height,width);
	dilation(m,height,width);
}

void bmpFile::open(unsigned char **m,int height,int width)
{
	dilation(m,height,width);
	erotion(m,height,width);
}

void bmpFile::extract(unsigned char **m,int height,int width)
{
	unsigned char** odata;
	odata = (unsigned char**)malloc(sizeof(unsigned char*)*bmpInfoHeader.biHeight);
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
	{
		odata[i]=(unsigned char*)malloc(nBytesPerLine);
		memcpy(odata[i],bmpFileData[i],nBytesPerLine);
	}
	dilation(m,height,width);
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		for(int j=0;j<bmpInfoHeader.biWidth;j++)
			if((bmpFileData[i][j] == 255) && (odata[i][j] == 0)) bmpFileData[i][j]=0;
			else bmpFileData[i][j]=255;
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		free(odata[i]);
	free(odata);
}

void bmpFile::HistogramLinearEqualizationWithLuminance()
{
#ifdef DEBUG
	cout << "HistogramLogEqualizationWithLuminance start.." << endl;
#endif
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	int min = 270, max=-1;
	memset(sumr,0,sizeof(sumr));
	memset(r,0,sizeof(r));
	if(bmpInfoHeader.biBitCount==24)
	{
		transferDataToYUV();
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[yuv[i][j*3]]++;
				if(yuv[i][j*3] < min) min = yuv[i][j*3];
				if(yuv[i][j*3] > max) max = yuv[i][j*3];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
			r[i] = int((sumr[i]*1.0)/tot*255);
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				yuv[i][j*3] = r[yuv[i][j*3]];
		transferYUVToData();
	} else {
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j]]++;
				if(bmpFileData[i][j] < min) min = bmpFileData[i][j];
				if(bmpFileData[i][j] > max) max = bmpFileData[i][j];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
			r[i] = int((sumr[i]*1.0)/tot*255);
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j] = r[bmpFileData[i][j]];
	}
#ifdef DEBUG
	cout << "HistogramLogEqualizationWithLuminance finished.." << endl;
#endif
}


void bmpFile::HistogramLinearEqualization()
{
#ifdef DEBUG
	cout << "HistogramLinearEqualization start.." << endl;
#endif
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	int max,min;
	if(bmpInfoHeader.biBitCount!=24) return;
	for(int k=0;k<3;k++)
	{
		max = -1;
		min = 256;
		memset(sumr,0,sizeof(sumr));
		memset(r,0,sizeof(r));
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j*3+k]]++;
				if(max < bmpFileData[i][j*3+k] ) max = bmpFileData[i][j*3+k];
				if(min > bmpFileData[i][j*3+k] ) min = bmpFileData[i][j*3+k];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
			r[i] = int((sumr[i]*1.0)/tot*255);
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*3+k] = r[bmpFileData[i][j*3+k]];
	}
#ifdef DEBUG
	cout << "HistogramLinearEqualization finished.." << endl;
#endif
}

void bmpFile::HistogramLogEqualizationWithLuminance(float min_pixel,float max_pixel,float base)
{
#ifdef DEBUG
	cout << "HistogramLogEqualizationWithLuminance start.." << endl;
#endif
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	int min = 270, max=-1;
	float b = (max_pixel-min_pixel)/(log(tot)/log(base));
	memset(sumr,0,sizeof(sumr));
	memset(r,0,sizeof(r));
	if(bmpInfoHeader.biBitCount==24)
	{
		transferDataToYUV();
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[yuv[i][j*3]]++;
				if(yuv[i][j*3] < min) min = yuv[i][j*3];
				if(yuv[i][j*3] > max) max = yuv[i][j*3];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
			r[i] = int(min_pixel+b*(log(sumr[i])/log(base))+0.5);
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				yuv[i][j*3] = r[yuv[i][j*3]];
		transferYUVToData();
	} else {
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j]]++;
				if(bmpFileData[i][j] < min) min = bmpFileData[i][j];
				if(bmpFileData[i][j] > max) max = bmpFileData[i][j];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
			r[i] = int(min_pixel+b*(log(sumr[i])/log(base))+0.5);
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j] = r[bmpFileData[i][j]];
	}
#ifdef DEBUG
	cout << "HistogramLogEqualizationWithLuminance finished.." << endl;
#endif
}

void bmpFile::HistogramLogEqualization(float min_pixel,float max_pixel,float base)
{
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	float b = (max_pixel-min_pixel)/(log(tot)/log(base));
	cout << " b=" << b <<endl;
	cout << " log 255 = " << log(255) << endl;
	cout << " log base = " << log(base) << endl;
	int max,min;
	if(bmpInfoHeader.biBitCount!=24) return;
	for(int k=0;k<3;k++)
	{
		max = -1;
		min = 256;
		memset(sumr,0,sizeof(sumr));
		memset(r,0,sizeof(r));
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j*3+k]]++;
				if(max < bmpFileData[i][j*3+k] ) max = bmpFileData[i][j*3+k];
				if(min > bmpFileData[i][j*3+k] ) min = bmpFileData[i][j*3+k];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
		{
			r[i] = int(min_pixel+b*(log(sumr[i])/log(base))+0.5);
			cout << i << "-> " << r[i] << endl;
		}
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*3+k] = r[bmpFileData[i][j*3+k]];
	}
}

void bmpFile::HistogramLogEqualization2(float a,float b,float base)
{
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	cout << " b=" << b <<endl;
	cout << " log 255 = " << log(255) << endl;
	cout << " log base = " << log(base) << endl;
	int max,min;
	if(bmpInfoHeader.biBitCount!=24) return;
	for(int k=0;k<3;k++)
	{
		max = -1;
		min = 256;
		memset(sumr,0,sizeof(sumr));
		memset(r,0,sizeof(r));
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j*3+k]]++;
				if(max < bmpFileData[i][j*3+k] ) max = bmpFileData[i][j*3+k];
				if(min > bmpFileData[i][j*3+k] ) min = bmpFileData[i][j*3+k];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
		{
			r[i] = int(a+b*(log(sumr[i]+1)/log(base)));
			if(r[i]>255) r[i] = 255;
			if(r[i]<0) r[i] = 0;
			cout << i << "-> " << r[i] << endl;
		}
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*3+k] = r[bmpFileData[i][j*3+k]];
	}
}

void bmpFile::HistogramLogEqualizationWithLuminance2(float a,float b,float base)
{
	int sumr[256];
	int r[256];
	int tot = bmpInfoHeader.biWidth*bmpInfoHeader.biHeight;
	int min = 270, max=-1;
	memset(sumr,0,sizeof(sumr));
	memset(r,0,sizeof(r));
	if(bmpInfoHeader.biBitCount==24)
	{
		transferDataToYUV();
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[yuv[i][j*3]]++;
				if(yuv[i][j*3] < min) min = yuv[i][j*3];
				if(yuv[i][j*3] > max) max = yuv[i][j*3];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
		{
			r[i] = int(a+b*(log(sumr[i]+1)/log(base)));
			if(r[i]>255) r[i] = 255;
			if(r[i]<0) r[i] = 0;
		}
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				yuv[i][j*3] = r[yuv[i][j*3]];
		transferYUVToData();
	} else {
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sumr[bmpFileData[i][j]]++;
				if(bmpFileData[i][j] < min) min = bmpFileData[i][j];
				if(bmpFileData[i][j] > max) max = bmpFileData[i][j];
			}
		for(int i=1;i<256;i++)
			sumr[i]+=sumr[i-1];
		for(int i=0;i<256;i++)
		{
			r[i] = int(a+b*(log(sumr[i]+1)/log(base)));
			if(r[i]>255) r[i] = 255;
			if(r[i]<0) r[i] = 0;
		}
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j] = r[bmpFileData[i][j]];
	}
}

void bmpFile::fill(int row,int col)
{

}

void bmpFile::mapResultMatrixWithCanvasSize(float canvasScaleCoef)
{
#ifdef DEBUG
	cout << "mapResultMatrix start.." << endl;
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "reversTransformMatrix = " << endl;
	reversTransformMatrix->printInfo("\t");
#endif
	int newWidth,newHeight;
	int oriHeight = bmpInfoHeader.biHeight;
	int oriWidth = bmpInfoHeader.biWidth;
	/*
	int minX,minY,maxX,maxY;

	Matrix* canvasParaMatrix = new Matrix(3,4);
	canvasParaMatrix->data[0][1] = bmpInfoHeader.biWidth-1;
	canvasParaMatrix->data[1][2] = bmpInfoHeader.biHeight-1;
	canvasParaMatrix->data[0][3] = bmpInfoHeader.biWidth-1;
	canvasParaMatrix->data[1][3] = bmpInfoHeader.biHeight-1;
	canvasParaMatrix->data[2][0] = canvasParaMatrix->data[2][1] = 1;
	canvasParaMatrix->data[2][2] = canvasParaMatrix->data[2][3] = 1;	

	canvasParaMatrix->printInfo("-------");

	canvasParaMatrix = Matrix::Multiplication(transformMatrix,canvasParaMatrix);

	canvasParaMatrix->printInfo("--------");
	Matrix* m = Matrix::Multiplication(new Matrix(*reversTransformMatrix),new Matrix(*canvasParaMatrix));

	m->printInfo("----------");
	maxX = minX = canvasParaMatrix->data[0][0];
	maxY = minY = canvasParaMatrix->data[1][0];
	for(int i=0;i<4;i++)
	{
		if(minX>canvasParaMatrix->data[0][i]) minX = canvasParaMatrix->data[0][i];
		if(maxX<canvasParaMatrix->data[0][i]) maxX = canvasParaMatrix->data[0][i];
		if(minY>canvasParaMatrix->data[1][i]) minY = canvasParaMatrix->data[1][i];
		if(maxY<canvasParaMatrix->data[1][i]) maxY = canvasParaMatrix->data[1][i];
	}

	cout << "-------------------------" << endl;
	cout << "MaxX = " << maxX << endl;
	cout << "MinX = " << minX << endl;
	cout << "MaxY = " << maxY << endl;
	cout << "minY = " << minY << endl;
	cout << "-------------------------" << endl;

	newWidth = int(maxX-minX+1);
	newHeight = int(maxY-minY+1);
	delete canvasParaMatrix;
	*/
	newWidth = bmpInfoHeader.biWidth*canvasScaleCoef;
	newHeight = bmpInfoHeader.biHeight*canvasScaleCoef;
	if(bmpInfoHeader.biBitCount == 24)
	{
#ifdef DEBUG
		cout << "bmpInfoHeader.biBitCount == 24" << endl;
		printInfo();
		cout << "ori width = " << bmpInfoHeader.biWidth << endl;
		cout << "ori height = " << bmpInfoHeader.biHeight << endl;
		cout << "ori nBytesPerLine = " << nBytesPerLine << endl;
#endif
		bmpInfoHeader.biWidth = newWidth;
		bmpInfoHeader.biHeight = newHeight;
		nBytesPerLine=(bmpInfoHeader.biWidth*3+3)/4*4;
		bmpInfoHeader.biSizeImage = nBytesPerLine*bmpInfoHeader.biHeight;
		bmpFileHeader.bfOffBits = 54;
		bmpFileHeader.bfSize = bmpFileHeader.bfOffBits + bmpInfoHeader.biSizeImage;
		unsigned char** data;
		data = (unsigned char**)malloc(bmpInfoHeader.biHeight*sizeof(unsigned char*));
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
		{
			data[i] = (unsigned char*)malloc(nBytesPerLine);
			memset(data[i],0,nBytesPerLine);
		}
		Matrix* resultPos = new Matrix(3,newHeight*newWidth);
		Matrix* originalPos;
		int row,col;
		for(int i=0;i<newHeight;i++)
			for(int j=0;j<newWidth;j++)
			{
				/*
				resultPos->data[0][i*newWidth+j] = j+minX;
				resultPos->data[1][i*newWidth+j] = i+minY;
				resultPos->data[2][i*newWidth+j] = 1;
				*/
				resultPos->data[0][i*newWidth+j] = j;
				resultPos->data[1][i*newWidth+j] = i;
				resultPos->data[2][i*newWidth+j] = 1;
			}
		originalPos = Matrix::Multiplication(reversTransformMatrix,resultPos);
		int rown,coln,a,b,c,d,m1,m2,m3,m4;
		float x,y;
		Matrix* m,*ans;
		m = new Matrix(4,5);
		for(int i=0;i<newHeight;i++)
			for(int j=0;j<newWidth;j++)
			{
				row = int(originalPos->data[1][i*newWidth+j]);
				col = int(originalPos->data[0][i*newWidth+j]);
				x = originalPos->data[0][i*newWidth+j]-col;
				y = originalPos->data[1][i*newWidth+j]-row;
				rown = row+1; coln = col+1;
				if(rown>=oriHeight) rown = row;
				if(coln>=oriWidth) coln = col;
				for(int k=0;k<3;k++)
					if(row<0||col<0||row>=oriHeight||col>=oriWidth) data[i][j*3+k] = 255;
					else 
					{
						/*
						a = bmpFileData[row][col*3+k];
						b = bmpFileData[rown][col*3+k];
						c = bmpFileData[row][coln*3+k];
						d = bmpFileData[rown][coln*3+k];
						m->data[0][0] = col; m->data[0][1] = row; m->data[0][2] = row*col; m->data[0][3] = 1; m->data[0][4] = a;
						m->data[0][0] = coln; m->data[0][1] = row; m->data[0][2] = row*coln; m->data[0][3] = 1;m->data[0][4] = c;
						m->data[0][0] = col; m->data[0][1] = rown; m->data[0][2] = rown*col; m->data[0][3] = 1;m->data[0][4] = b;
						m->data[0][0] = coln; m->data[0][1] = rown; m->data[0][2] = rown*coln; m->data[0][3] = 1;m->data[0][4] = d;
						ans = m->SolveEquation();
						if(ans!=NULL)
						{
							a = ans->data[0][0];
							b = ans->data[1][0];
							c = ans->data[2][0];
							d = ans->data[3][0];
							d = a*x+b*y+c*x*y+d;
							if(d>255)
								data[i][j*3+k] = 255;
							else data[i][j*3+k] = int(d);
						} else data[i][j*3+k] = bmpFileData[row][col*3+k];
						*/
						d = bmpFileData[row][col*3+k];
						a = bmpFileData[row][coln*3+k]-d;
						b = bmpFileData[rown][col*3+k]-d;
						c = bmpFileData[rown][coln*3+k]-a-b-d;
						d = a*x+b*y+c*x*y+d;
						if(d>255)
							data[i][j*3+k] = 255;
						else if(d<0)
						   	data[i][j*3+k] = 0;
						else data[i][j*3+k] = int(d);
					}
			}
		reversTransformMatrix = Matrix::UnitMatrix(3);
		transformMatrix = Matrix::UnitMatrix(3);
		for(int i=0;i<oriHeight;i++)
			free(bmpFileData[i]);
		free(bmpFileData);
		bmpFileData = data;
	}
#ifdef DEBUG
	cout << "after width = " << bmpInfoHeader.biWidth << endl;
	cout << "after height = " << bmpInfoHeader.biHeight << endl;
	cout << "after nBytesPerLine = " << nBytesPerLine << endl;
	printInfo();
	cout << "mapResultMatrix finished.." << endl;
#endif
}

void bmpFile::mirrorX(int x)
{
#ifdef DEBUG
	cout << "mirrorX start.." << endl;
#endif
	transformMatrix = Matrix::MirrorX(transformMatrix,x);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::MirrorX(m,x);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "mirrorX finished.." << endl;
#endif
}

void bmpFile::mirrorY(int y)
{
#ifdef DEBUG
	cout << "mirrorY start.." << endl;
#endif
	transformMatrix = Matrix::MirrorY(transformMatrix,y);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::MirrorY(m,y);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "mirrorY finished.." << endl;
#endif
}

void bmpFile::translate(int x,int y)
{
#ifdef DEBUG
	cout << "translate start.." << endl;
#endif
	transformMatrix = Matrix::Translate(transformMatrix,x,y);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::Translate(m,-x,-y);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "translate finished.." << endl;
#endif
}

void bmpFile::scale(float x,float y)
{
#ifdef DEBUG
	cout << "scale start.." << endl;
#endif
	transformMatrix = Matrix::Scale(transformMatrix,x,y);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::Scale(m,1/x,1/y);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "scale finished.." << endl;
#endif
}

void bmpFile::rotate(float degree) 
{
#ifdef DEBUG
	cout << "rotate start.." << endl;
#endif
	transformMatrix = Matrix::Rotate(transformMatrix,degree);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::Rotate(m,-degree);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "Rotate finished.." << endl;
#endif
}

void bmpFile::shear(float x,float y) 
{
#ifdef DEBUG
	cout << "shear start.." << endl;
#endif
	transformMatrix = Matrix::Shear(transformMatrix,x,y);

	Matrix* m =Matrix::UnitMatrix(3);
	m = Matrix::Shear(m,-x,-y);
	reversTransformMatrix = Matrix::Multiplication(reversTransformMatrix,m);
#ifdef DEBUG
	cout << "transformMatrix = " << endl;
	transformMatrix->printInfo("\t");
	cout << "Rotate finished.." << endl;
#endif
}

unsigned char bmpFile::clip(int k)
{
	if(k<0) k = 0;
	if(k>255) k = 255;
	return (unsigned char)k;
}

unsigned char bmpFile::clip(float k)
{
	if(k<0) k = 0;
	if(k>255) k = 255;
	return (unsigned char)k;
}

float bmpFile::Convolution(Matrix* filter,int r,int c,int times,int k)
{
	float res = 0;
	int startRow = r-filter->Row/2;
	int startCol = c-filter->Col/2;
	int row,col;
	for(int i=0;i<filter->Row;i++)
	{
		row = i+startRow;
		if(row>=bmpInfoHeader.biHeight) break;
		if(row<0) continue;
		for(int j=0;j<filter->Col;j++)
		{
			col = startCol+j;
			if(col>=bmpInfoHeader.biWidth) break;
			if(col<0) continue;
			res+= filter->data[i][j]*bmpFileData[row][col*times+k];
		}
	}
	return res;
}

void bmpFile::MeanFilter(int length)
{
	float sum = length*length;
	Matrix* filter = Matrix::ones(length);
	Matrix m[3] = {
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth)
	};
	int times = 1;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	for(int k=0;k<times;k++)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				m[k][i][j] = Convolution(filter,i,j,times,k)/sum;	
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = fabs(bmpFileData[i][j*times+k]-(unsigned char)m[k][i][j])*3;
	}
	exportToFile("!!.bmp");
	for(int k=0;k<times;k++)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)m[k][i][j];
	}
	delete filter;
}

void bmpFile::GaussianFilter(float sigma)
{
	Matrix* filter = Matrix::GaussianFilter(sigma);
	Matrix m[3] = {
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth)
	};
	float sum = 0;
	int times = 1;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	for(int i=0;i<filter->Row;i++)
		for(int j=0;j<filter->Col;j++)
			sum+= filter->data[i][j];
	for(int k=0;k<times;k++)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				m[k][i][j] = Convolution(filter,i,j,times,k)/sum;	 
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = fabs(bmpFileData[i][j*times+k]-(unsigned char)m[k][i][j])*3;
	}
//	exportToFile("!!.bmp");
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)m[k][i][j];
	delete filter;
}

void bmpFile::generateBilateralFilter(Matrix* inputFilter,Matrix* resultFilter,float sigmar,int r,int c,int times,int k)
{
	int startRow = r-inputFilter->Row/2;
	int startCol = c-inputFilter->Col/2;
	int row,col;
	float intensity_dis;
	for(int i=0;i<inputFilter->Row;i++)
	{
		row = i+startRow;
		if(row>=bmpInfoHeader.biHeight) break;
		if(row<0) continue;
		for(int j=0;j<inputFilter->Col;j++)
		{
			col = startCol+j;
			if(col>=bmpInfoHeader.biWidth) break;
			if(col<0) continue;
			intensity_dis = fabs(bmpFileData[row][col*times+k]-bmpFileData[r][c*times+k]);
			resultFilter->data[i][j] = Matrix::GaussianFunction(sigmar,intensity_dis)*inputFilter->data[i][j];
		}
	}
}

void bmpFile::generateColorfulBilateralFilter(Matrix* inputFilter,Matrix* resultFilter,float sigmar,int r,int c)
{
	int startRow = r-inputFilter->Row/2;
	int startCol = c-inputFilter->Col/2;
	int row,col;
	float intensity_dis;
	int times = 3;
	for(int i=0;i<inputFilter->Row;i++)
	{
		row = i+startRow;
		if(row>=bmpInfoHeader.biHeight) break;
		if(row<0) continue;
		for(int j=0;j<inputFilter->Col;j++)
		{
			col = startCol+j;
			if(col>=bmpInfoHeader.biWidth) break;
			if(col<0) continue;
			intensity_dis = 0;
			for(int k=0;k<times;k++)
				intensity_dis += sqr(bmpFileData[row][col*times+k]-bmpFileData[r][c*times+k]);
			intensity_dis = sqrt(intensity_dis);
			resultFilter->data[i][j] = Matrix::GaussianFunction(sigmar,intensity_dis)*inputFilter->data[i][j];
		}
	}
}

void bmpFile::BilateralFilter(float sigmas,float sigmar)
{
	Matrix* filter = Matrix::GaussianFilter(sigmas);
	Matrix* GaussianFilter = Matrix::GaussianFilter(sigmas);
	Matrix m[3] = {
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth)
	};

	float sum[3];
	float min[3]={1000,1000,1000},max[3] = {-1000,-1000,-1000};
	int times;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	else times = 1;
	for(int k=0;k<times;k++)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				sum[k] = 0;
				generateBilateralFilter(GaussianFilter,filter,sigmar,i,j,times,k);
				for(int ii = 0;ii<filter->Row;ii++)
					for(int jj=0;jj<filter->Col;jj++)
						sum[k]+=filter->data[ii][jj];
				m[k][i][j] = Convolution(filter,i,j,times,k);
				if(sum[k]!=0) m[k][i][j] = m[k][i][j]/sum[k];
			}
		/*
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				if(min[k] > bmpFileData[i][j*times+k]-m[k][i][j]) min[k] =(bmpFileData[i][j*times+k]-m[k][i][j]);
				if(max[k] < bmpFileData[i][j*times+k]-m[k][i][j]) max[k] =(bmpFileData[i][j*times+k]-m[k][i][j]);
			}
		float gap = max - min;
		*/
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)clip(bmpFileData[i][j*times+k]-m[k][i][j]);
	}
//	exportToFile("!!.bmp");
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)m[k][i][j];
}

void bmpFile::ColorfulBilateralFilter(float sigmas,float sigmar)
{
	Matrix* filter = Matrix::GaussianFilter(sigmas);
	Matrix* GaussianFilter = Matrix::GaussianFilter(sigmas);
	Matrix m[3] = {
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth)
	};

	float sum;
	int times;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	else return;
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		for(int j=0;j<bmpInfoHeader.biWidth;j++)
		{
			sum = 0;
			generateColorfulBilateralFilter(GaussianFilter,filter,sigmar,i,j);
			for(int ii = 0;ii<filter->Row;ii++)
				for(int jj=0;jj<filter->Col;jj++)
					sum+=filter->data[ii][jj];
			for(int k=0;k<times;k++)
			{
				m[k][i][j] = Convolution(filter,i,j,times,k);
				if(sum!=0) m[k][i][j] = m[k][i][j]/sum;
			}
		}
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = clip((bmpFileData[i][j*times+k]-m[k][i][j])*10);
//	exportToFile("!!.bmp");
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)m[k][i][j];
}

float bmpFile::getVariance(Matrix* inputFilter,int r,int c,int times)
{
	int startRow = r-inputFilter->Row/2;
	int startCol = c-inputFilter->Col/2;
	int row,col;
	float sum[3] = {0,0,0};
	float number = 0;
	for(int i=0;i<inputFilter->Row;i++)
	{
		row = i+startRow;
		if(row>=bmpInfoHeader.biHeight) break;
		if(row<0) continue;
		for(int j=0;j<inputFilter->Col;j++)
		{
			col = startCol+j;
			if(col>=bmpInfoHeader.biWidth) break;
			if(col<0) continue;
			number++;
			for(int k=0;k<times;k++)
				sum[k]+= bmpFileData[row][col*times+k];
		}
	}
	float average [3];
	for(int k=0;k<times;k++)
		average[k] = sum[k]/number;
	float res = 0,tmp;
	for(int i=0;i<inputFilter->Row;i++)
	{
		row = i+startRow;
		if(row>=bmpInfoHeader.biHeight) break;
		if(row<0) continue;
		for(int j=0;j<inputFilter->Col;j++)
		{
			col = startCol+j;
			if(col>=bmpInfoHeader.biWidth) break;
			if(col<0) continue;
			tmp = 0;
			for(int k=0;k<times;k++)
				tmp += sqr(bmpFileData[row][col*times+k]-average[k]);
			res += sqrt(tmp);
		}
	}
//	return sqrt(res)/(number-1);
	return res/number;
}

void bmpFile::AutoBilateralFilter(float sigmas)
{
	float sigmar = 0;
	Matrix* filter = Matrix::GaussianFilter(sigmas);
	Matrix* GaussianFilter = Matrix::GaussianFilter(sigmas);
	Matrix m[3] = {
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth),
		Matrix(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth)
	};

	float sum;
	int times;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	else times = 1;
	for(int i=0;i<bmpInfoHeader.biHeight;i++)
		for(int j=0;j<bmpInfoHeader.biWidth;j++)
		{
			sum = 0;
			sigmar = getVariance(GaussianFilter,i,j,times);
			generateColorfulBilateralFilter(GaussianFilter,filter,sigmar,i,j);
			for(int ii = 0;ii<filter->Row;ii++)
				for(int jj=0;jj<filter->Col;jj++)
					sum+=filter->data[ii][jj];
			for(int k=0;k<times;k++)
			{
				m[k][i][j] = Convolution(filter,i,j,times,k);
				if(sum!=0) m[k][i][j] = m[k][i][j]/sum;
			}
		}
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = clip((bmpFileData[i][j*times+k]-m[k][i][j])*10);
//	exportToFile("!!.bmp");
	for(int k=0;k<times;k++)
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = (unsigned char)m[k][i][j];
}

void bmpFile::LaplacianEnhancement()
{
	int length = 3;
	Matrix* filter = Matrix::Laplacian(length);
	Matrix m(bmpInfoHeader.biHeight,bmpInfoHeader.biWidth);
	int times = 1;
	float min = 10000,max = -10000;
	float min_ori = 10000,max_ori = -10000;
	float gap;
	float range;
	if(bmpInfoHeader.biBitCount == 24) times = 3;
	for(int k=0;k<times;k++)
	{
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
			{
				m[i][j] = Convolution(filter,i,j,times,k);
				if(min>m[i][j]+bmpFileData[i][j*times+k]) min = m[i][j]+bmpFileData[i][j*times+k];
				if(max<m[i][j]+bmpFileData[i][j*times+k]) max = m[i][j]+bmpFileData[i][j*times+k];
				if(min_ori>bmpFileData[i][j*times+k]) min_ori = bmpFileData[i][j*times+k];
				if(max_ori<bmpFileData[i][j*times+k]) max_ori = bmpFileData[i][j*times+k];
			}
		gap = max - min;
//		range = (max_ori-min_ori)/2;
		range = 255;
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
//				m[i][j] = (m[i][j]-min)/gap*range*2-range;
				if(m[i][j]<0) m[i][j] = -m[i][j]/min*range;
				else m[i][j] = m[i][j]/max*range;
		for(int i=0;i<bmpInfoHeader.biHeight;i++)
			for(int j=0;j<bmpInfoHeader.biWidth;j++)
				bmpFileData[i][j*times+k] = clip(int(m[i][j]+bmpFileData[i][j*times+k]));
	}
}
