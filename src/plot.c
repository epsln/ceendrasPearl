#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "include/plot.h"

int checkBoundaries(int x, int y, image_t *img){
	if (x >= 0 && y >= 0 && x < img->w && y < img->h)
		return 1;
	else
		return 0;
}

void plotLineLow(int x0,int y0, int x1,int y1, image_t* img){
	int dx = x1 - x0;
	int dy = y1 - y0;
	int yi = 1;
	if (dy < 0){
		yi = -1;
		dy = -dy;
	}
	int D = (2 * dy) - dx;
	int y = y0;

	for (int x = x0; x < x1; x++){
		point(x,y, img);
		if (D > 0){
			y = y + yi;
			D = D + (2 * (dy - dx));
		}
		else
			D = D + 2*dy;
	}
}

void plotLineHigh(int x0,int y0, int x1,int y1, image_t* img){
	int dx = x1 - x0;
	int dy = y1 - y0;
	int xi = 1;
	if (dx < 0){
		xi = -1;
		dx = -dx; 
	}
	int D = (2 * dx) - dy;
	int x = x0;

	for (int y = y0; y < y1; y++){	
		point(x,y, img);
		if (D > 0){
			x = x + xi;
			D = D + (2 * (dx - dy));
		}
		else
			D = D + 2*dx;
	}
}

void point(int x, int y, image_t* img){
	if (checkBoundaries(x, y, img) == 0) return;
	else{
		img->pointArr[x][y] = 1;
	}
}

void line(int x0,int y0, int x1, int y1, image_t* img){
	if (img->line == 0) return;

	if (checkBoundaries(x0, y0, img) == 0 || checkBoundaries(x1, y1, img) == 0 ) return;

	if (abs(y1 - y0) < abs(x1 - x0)){
		if (x0 > x1)
			plotLineLow(x1, y1, x0, y0, img);
		else
			plotLineLow(x0, y0, x1, y1, img);
	}
	else{
		if (y0 > y1)
			plotLineHigh(x1, y1, x0, y0, img);
		else
			plotLineHigh(x0, y0, x1, y1, img);
	}
}

void antialiasing(image_t* img, float*** output){
	const int w0 = img->w;
	const int h0 = img->h;

	const int antPow = img->antialiasingPow;
	for (int i = 0; i < w0; i++){
		for (int j = 0; j < h0; j++){
			output[i/antPow][j/antPow][0] += img->pointArr[i][j]; 	
			output[i/antPow][j/antPow][1] += img->pointArr[i][j]; 	
			output[i/antPow][j/antPow][2] += img->pointArr[i][j]; 	

			output[i/antPow][j/antPow][0] /= 2;
			output[i/antPow][j/antPow][1] /= 2;
			output[i/antPow][j/antPow][2] /= 2;
		}
	}
}

void saveArrayAsBMP(image_t *img){
	int w = img->w/img->antialiasingPow;
	int h = img->h/img->antialiasingPow;

	FILE *f;
	unsigned char *imgOut = NULL;
	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int
	//Allocate image array
	float *** imgArr = NULL;
	imgArr = (float***)malloc(w*sizeof(float**));
	for (int i = 0; i < w; i++) {
		imgArr[i] = (float **) malloc(h*sizeof(float *));
		for (int j = 0; j < h; j++){
			imgArr[i][j] = (float *) malloc(3 *sizeof(float));
			imgArr[i][j][0] = 0;
			imgArr[i][j][1] = 0;
			imgArr[i][j][2] = 0;
		}
	}

	if (imgArr == NULL){
		printf("Could not allocate memory for the image array !\nExiting...\n");
		exit(-1);
	}

	antialiasing(img, imgArr);	

	imgOut = (unsigned char *)malloc(3*w*h);

	if (imgOut == NULL){
		printf("ERROR: Failed to allocate memory for the image !!\nExiting...\n");
		exit(1);
	}

	memset(imgOut,0,3*w*h);

	int r,g,b,x,y;
	for(int i=0; i<w; i++)
	{
		for(int j=0; j<h; j++)
		{
			x=i; y=(h-1)-j;
			r = imgArr[i][j][0]*255;
			g = imgArr[i][j][1]*255;
			b = imgArr[i][j][2]*255;
			if (r > 255) r=255;
			if (g > 255) g=255;
			if (b > 255) b=255;
			imgOut[(x+y*w)*3+2] = (unsigned char)(r);
			imgOut[(x+y*w)*3+1] = (unsigned char)(g);
			imgOut[(x+y*w)*3+0] = (unsigned char)(b);
		}
	}

	unsigned char bmpfileheader[14] = {'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0};
	unsigned char bmpinfoheader[40] = {40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0};
	unsigned char bmppad[3] = {0,0,0};

	bmpfileheader[ 2] = (unsigned char)(filesize    );
	bmpfileheader[ 3] = (unsigned char)(filesize>> 8);
	bmpfileheader[ 4] = (unsigned char)(filesize>>16);
	bmpfileheader[ 5] = (unsigned char)(filesize>>24);

	bmpinfoheader[ 4] = (unsigned char)(       w    );
	bmpinfoheader[ 5] = (unsigned char)(       w>> 8);
	bmpinfoheader[ 6] = (unsigned char)(       w>>16);
	bmpinfoheader[ 7] = (unsigned char)(       w>>24);
	bmpinfoheader[ 8] = (unsigned char)(       h    );
	bmpinfoheader[ 9] = (unsigned char)(       h>> 8);
	bmpinfoheader[10] = (unsigned char)(       h>>16);
	bmpinfoheader[11] = (unsigned char)(       h>>24);

	f = NULL;
	f = fopen(img->filename,"wb");
	if (f == NULL){
		printf("Could not open image... Does the output folder exists ?\nExiting...\n");
		exit(2);
	}	
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(int i=0; i<h; i++)
	{
		fwrite(imgOut+(w*(h-i-1)*3),3,w,f);
		fwrite(bmppad,1,(4-(w*3)%4)%4,f);
	}

	//Free the memory 
	for (int i = 0; i < w; i++) {
		for (int j = 0; j < h; j++){
			free(imgArr[i][j]);
		}
		free(imgArr[i]);
	}
	free(imgArr);

	free(imgOut);
	fclose(f);
}


