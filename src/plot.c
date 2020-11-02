#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


typedef struct image{
	int width;
	int height;

	int** hitArray; //array which will store where the iterated points end up 

	float*** imgArray; //Actual image; can be downscaled for antialiasing
}image;

int checkBoundaries(int x, int y, int maxW, int maxH){
	if (x >= 0 && y >= 0 && x < maxW && y < maxH)
		return 1;
	else
		return 0;
}

void plotLineLow(int x0,int y0, int x1,int y1, float*** imgArr){
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
		imgArr[x][y][0] = 1;
		imgArr[x][y][1] = 1;
		imgArr[x][y][2] = 1;
		if (D > 0){
			y = y + yi;
			D = D + (2 * (dy - dx));
		}
		else
			D = D + 2*dy;
	}
}

void plotLineHigh(int x0,int y0, int x1,int y1, float*** imgArr){
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
		imgArr[x][y][0] = 1.;
		imgArr[x][y][1] = 1.;
		imgArr[x][y][2] = 1.;
		if (D > 0){
			x = x + xi;
			D = D + (2 * (dx - dy));
		}
		else
			D = D + 2*dx;
	}
}

void point(int x, int y, float*** imgArr, int w, int h){
	if (checkBoundaries(x, y, w, h) == 0) return;
	else{
		imgArr[x][y][0] = 1;
		imgArr[x][y][1] = 1;
		imgArr[x][y][2] = 1;
	}
}

void line(int x0,int y0, int x1,int y1, float*** imgArr, int LINE, int w, int h ){
	if (LINE == 0) return;

	if (checkBoundaries(x0, y0, w, h) == 0 || checkBoundaries(x1, y1, w, h) == 0 ) return;

	if (abs(y1 - y0) < abs(x1 - x0)){
		if (x0 > x1)
			plotLineLow(x1, y1, x0, y0, imgArr);
		else
			plotLineLow(x0, y0, x1, y1, imgArr);
	}
	else{
		if (y0 > y1)
			plotLineHigh(x1, y1, x0, y0, imgArr);
		else
			plotLineHigh(x0, y0, x1, y1, imgArr);
	}
}

void antialiasing(float*** imgArr, double* PARAMS, float*** output, int power){
	int w0 = PARAMS[3];
	int h0 = PARAMS[4];


	for (int i = 0; i < w0; i++){
		for (int j = 0; j < h0; j++){
			output[i/power][j/power][0] += imgArr[i][j][0]; 	
			output[i/power][j/power][1] += imgArr[i][j][1]; 	
			output[i/power][j/power][2] += imgArr[i][j][2]; 	
		}
	}
	for (int i = 0; i < w0/power; i++){
		for (int j = 0; j < h0/power; j++){

			output[i][j][0] /= pow(2, power);
			output[i][j][1] /= pow(2, power);
			output[i][j][2] /= pow(2, power);
		}
	}
}

void saveArrayAsBMP(float*** input, char* filename, double* PARAMS){
	int downsamplePower = 4; //Downsampling 4 times
	int w = PARAMS[3]/downsamplePower;
	int h = PARAMS[4]/downsamplePower;

	FILE *f;
	unsigned char *img = NULL;
	int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int

	//Allocate image array
	float *** imgArr = NULL;
	imgArr = (float***)malloc(w*sizeof(float**));
	for (int i = 0; i< w; i++) {
		imgArr[i] = (float **) malloc(h*sizeof(float *));
		for (int j = 0; j < h; j++) 
			imgArr[i][j] = (float *) malloc(3 *sizeof(float));
	}

	//Zero all elements
	for (int i = 0; i< w; i++) {
		for (int j = 0; j < h; j++){
			imgArr[i][j][0] = 0;
			imgArr[i][j][1] = 0;
			imgArr[i][j][2] = 0;
		}
	}

	if (imgArr == NULL){
		printf("Could not allocate memory for the image array !\nExiting...\n");
		exit(-1);
	}

	antialiasing(input, PARAMS, imgArr, downsamplePower);	

	img = (unsigned char *)malloc(3*w*h);
	memset(img,0,3*w*h);
	if (img == NULL){
		printf("ERROR: Failed to allocate memory for the image !!\nExiting...\n");
		exit(1);
	}
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
			img[(x+y*w)*3+2] = (unsigned char)(r);
			img[(x+y*w)*3+1] = (unsigned char)(g);
			img[(x+y*w)*3+0] = (unsigned char)(b);
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
	f = fopen(filename,"wb");
	if (f == NULL){
		printf("Could not open image... Does the output folder exists ?\nExiting...\n");
		exit(2);
	}	
	fwrite(bmpfileheader,1,14,f);
	fwrite(bmpinfoheader,1,40,f);
	for(int i=0; i<h; i++)
	{
		fwrite(img+(w*(h-i-1)*3),3,w,f);
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

	free(img);
	fclose(f);
}


