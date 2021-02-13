#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "include/plot.h"
#include "include/complexMath.h"

int checkBoundaries(int x, int y, image_t *img){
    if (img->bitwise){
        if (x >= 0 && y >= 0 && x < img->w && y < img->h &&  x*img->h + y > 0) //Careful, might overflow !
            return 1;
        else
            return 0;

    }
    if (x >= 0 && y >= 0 && x*img->h + y < img->w * img->h &&  x*img->h + y > 0) //Careful, might overflow !
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

unsigned long long output(unsigned long long n)
{
    unsigned long long m = n ? output(n / 2) : 0;
    printf("%d", (int)(n % 2));
    return m;
}

void point(int x, int y, image_t* img){
    if (checkBoundaries(x, y, img) == 0) return;
    if (img->bitwise == 1){
        img->bitArray[(int)fmax(0.0, (int)ceil(x/63.0) - 1) * img->h + (int)y] |= 1ULL << (int)(63 - x % 64) ;
    }
    else{
        img->pointArr[x*img->h + y] = 1;
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


void antialiasing(image_t* img, unsigned char* outputImg){
    const int w0 = img->w;
    const int h0 = img->h;

    const int antPow = img->antialiasingPow;
    const int minPixelValue = 255/(antPow * 2);
    //What we do here:
    //Copy each bit of the bitArray into a float at its corresponding index i,j
    //We use some bit shuffling to select the correct bit, starting from the end (we are storing bits in a big endian manner) 
    //we use i as the index for this bit. it runs from 63 to 0
    //Each time we add a bit, we also divide it by two to obtain some kind of mean

    if (img->bitwise == 1){
        for (int i = 0; i < h0; i++){
            for (int j = 0; j < w0; j++){
                int res = minPixelValue * ((img->bitArray[(int)fmax(0, ceil(j/63.0) - 1) * img->h + i] & (1ULL << (63 - j % 64))) >> (63 - j % 64));
                if (res == 0)
                    continue;
                outputImg[(i/antPow * w0/antPow + j/antPow) * 3 + 0] += res;
                outputImg[(i/antPow * w0/antPow + j/antPow) * 3 + 1] += res; 
                outputImg[(i/antPow * w0/antPow + j/antPow) * 3 + 2] += res;
            }
        }
        //Zero bit array after reading
	//Bugged !
        memset(img->bitArray, 0, (ceil(img->w/64.0) + 2)*img->h *(sizeof(long long int)));
    }
    //Classical method, just add up all the floats and then divide
    else{
        for (int i = 0; i < h0; i++) {
            for (int j = 0; j < w0; j++){
                int res = minPixelValue * img->pointArr[j * img->h + i];
                outputImg[(j/antPow + i/antPow * w0/antPow ) * 3 + 0] += res; 	
                outputImg[(j/antPow + i/antPow * w0/antPow ) * 3 + 1] += res; 	
                outputImg[(j/antPow + i/antPow * w0/antPow ) * 3 + 2] += res; 	
            }
        }
        memset(img->pointArr, 0, (img->w*img->h) * (sizeof *img->pointArr));
    }
}

void makeFilename(int numImg, image_t *img){
	char prefix[100] = "out/img_";
	char numStr[6];  

		sprintf(numStr, "%d", numImg);
		strcat(prefix, numStr);
		strcat(prefix, ".bmp\0");
		strcpy(img->filename, prefix);
		memset(prefix, 0, sizeof(prefix));
		strcpy(prefix, "out/img_");


}


void saveArrayAsBMP(image_t *img){
    int w = img->w/img->antialiasingPow;
    int h = img->h/img->antialiasingPow;

    FILE *f;
    int filesize = 54 + 3*w*h;  //w is your image width, h is image height, both int
    //Allocate image array
    float * imgArr = NULL;
    unsigned char *imgOut = NULL;

    imgArr = (float*)calloc(w * h * 3, sizeof(float));
    imgOut = (unsigned char *)calloc(3*w*h, sizeof(unsigned char));

    if (imgArr == NULL || imgOut == NULL){
        printf("Could not allocate memory for the image array !\nExiting...\n");
        exit(1);
    }

    antialiasing(img, imgOut);	

    //memset(img->pointArr,0,img->w*img->h);

    //int r,g,b,x,y;
    //for(int i=0; i<w; i++)
    //{
    //	for(int j=0; j<h; j++)
    //	{
    //		x=i; y=(h-1)-j;
    //		r = imgArr[(i + j * w) * 3 + 0]*255;
    //		g = imgArr[(i + j * w) * 3 + 1]*255;
    //		b = imgArr[(i + j * w) * 3 + 2]*255;
    //		if (r > 255) r=255;
    //		if (g > 255) g=255;
    //		if (b > 255) b=255;
    //		imgOut[(x+y*w)*3+2] = (unsigned char)(r);
    //		imgOut[(x+y*w)*3+1] = (unsigned char)(g);
    //		imgOut[(x+y*w)*3+0] = (unsigned char)(b);
    //	}
    //}

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

    //free(imgArr);

    free(imgOut);
    fclose(f);
}


