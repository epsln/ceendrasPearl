#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <complex.h>

#include "include/plot.h"
#include "include/complexMath.h"
#include "include/arraysOps.h"
#include "include/treeExploration.h"
#include "include/debugTools.h"

#define SIZEARR 1000
#define HEIGHT 1080 * 4
#define WIDTH  1920 * 4
#define BOUNDS 2 
#define EPSI  0.01 
#define LINE 0 
#define DEBUG 0



int main(){
        time_t t;
	srand((unsigned) time(&t));

	double complex ta = 0.;
	double complex tb = 0.;
	double complex taBeg = 0.;
	double complex tbBeg = 0.;
	double complex taEnd=  0.;
	double complex tbEnd = 0.;
	double complex tab = 0.;
	int numIm = 0;
	float theta = 0;
	taBeg = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
	tbBeg = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
	taEnd = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
	tbEnd = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
	printf("taBeg: %lf + %lf\n", creal(taBeg), cimag(taBeg));
	printf("taEnd: %lf + %lf\n", creal(taEnd), cimag(taEnd));
	printf("tbBeg: %lf + %lf\n", creal(tbBeg), cimag(tbBeg));
	printf("tbEnd: %lf + %lf\n", creal(tbEnd), cimag(tbEnd));
	srand(time(NULL));

	double PARAMS[10];
	//TODO: Use a config file for this ;)
	PARAMS[0] = 10; //Maximum depth
	PARAMS[1] = EPSI; //epsilon
	PARAMS[2] = BOUNDS; //bounds
	PARAMS[3] = WIDTH; //x resolution
	PARAMS[4] = HEIGHT; //y resolution
	PARAMS[5] = LINE; //Draw line (1 = yes)
	PARAMS[6] = DEBUG; //Debug mode (1 = yes)

	float *** imgArr = NULL;
	imgArr = (float***)malloc(WIDTH*sizeof(float**));
	for (int i = 0; i< WIDTH; i++) {
		imgArr[i] = (float **) malloc(HEIGHT*sizeof(float *));
		for (int j = 0; j < HEIGHT; j++) 
			imgArr[i][j] = (float *) malloc(3 *sizeof(float));
	}

	if (imgArr == NULL){
		printf("Could not allocate memory for the image array !\nExiting...\n");
		exit(-1);
	}
	while(1){

		for (int i = 0; i< WIDTH; i++) {
			for (int j = 0; j < HEIGHT; j++) {
				imgArr[i][j][0] = 0;
				imgArr[i][j][1] = 0;
				imgArr[i][j][2] = 0;
			}
		}

		ta = InOutQuadComplex((float)numIm, taBeg, taEnd, 300); 
		tb = InOutQuadComplex((float)numIm, tbBeg, tbEnd, 300);
		computeDepthFirst(PARAMS, ta, tb, tab, imgArr, numIm);
		numIm++;
		//	LEVMAX++;
		if (numIm > 300 ) return(1);
		if (theta > 3.1415928 + 3.1415928/(30*10) ) return(1);
		theta += 3.1415928/(30*10);
	}
	return 1;
}
