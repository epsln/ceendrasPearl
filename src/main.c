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
#define EPSI  0.005 
#define LINE 0 
#define DEBUG 0



int main(){
	time_t t;
	srand((unsigned) time(&t));

	double complex ta = 0.;
	double complex tb = 0.;
	double complex taBeg = 0.;
	double complex tbBeg = 0.;
	double complex taEnd =  0.;
	double complex tbEnd = 0.;
	double complex tab = 0.;
	double complex taInit = 0.;
	double complex tbInit = 0.;
	int numIm = 0;
	float theta = 0;
	taBeg = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
	tbBeg = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
	taEnd = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
	tbEnd = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
	taInit = taBeg;
	tbInit = tbBeg;
	printf("taBeg: %lf + %lf\n", creal(taBeg), cimag(taBeg));
	printf("taEnd: %lf + %lf\n", creal(taEnd), cimag(taEnd));
	printf("tbBeg: %lf + %lf\n", creal(tbBeg), cimag(tbBeg));
	printf("tbEnd: %lf + %lf\n\n", creal(tbEnd), cimag(tbEnd));
	srand(time(NULL));

	double PARAMS[10];
	//TODO: Use a config file for this ;)
	PARAMS[0] = 14; //Maximum depth
	PARAMS[1] = EPSI; //epsilon
	PARAMS[2] = BOUNDS; //bounds
	PARAMS[3] = WIDTH; //x resolution
	PARAMS[4] = HEIGHT; //y resolution
	PARAMS[5] = LINE; //Draw line (1 = yes)
	PARAMS[6] = DEBUG; //Debug mode (1 = yes)

	int fps = 30;
	int duration = 6;
	int lengthAnim = 30;
	float *** imgArr = NULL;
	imgArr = (float***)malloc(WIDTH*sizeof(float**));//TODO:Replace by calloc
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

		//Here, we interpolate between two values of ta,tb using an easing function
		//The easing function takes a starting value and the value that needs to be added
		//To get the value that needs to be added we extract the distance between the two traces using copysign
		//And multiply by minus one to add. We need to do that for the real and complex part so we get this loooong line :)
		ta = InOutQuadComplex((float)(numIm%(fps*duration)), taBeg, -copysign(creal(taBeg- taEnd), creal(taBeg- taEnd)) + I*-copysign(cimag(taBeg- taEnd), cimag(taBeg- taEnd)), (float)fps * duration); 
		tb = InOutQuadComplex((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		printf("%lf + i %lf\n", creal(ta), cimag(ta));
		printf("%lf + i %lf\n\n", creal(tb), cimag(tb));

		ta = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
		tb = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
		computeDepthFirst(PARAMS, ta, tb, tab, imgArr, numIm);
		numIm++;
		if (numIm % (fps * duration) == 0 ){
			taBeg = taEnd;
			tbBeg = tbEnd;
			if (numIm >= fps * lengthAnim - fps * duration){//loop by ending up at the begining traces
				printf("looping !\n");
				taBeg = taEnd;
				tbBeg = tbEnd;
				taEnd = taInit;
				tbEnd = tbInit;
			}
			else{
				taEnd = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
				tbEnd = randomComplex(-3 + 1.5 * I, 3 + 1.5 * I);
			}
			printf("taBeg: %lf + %lf\n", creal(taBeg), cimag(taBeg));
			printf("taEnd: %lf + %lf\n", creal(taEnd), cimag(taEnd));
			printf("tbBeg: %lf + %lf\n", creal(tbBeg), cimag(tbBeg));
			printf("tbEnd: %lf + %lf\n\n", creal(tbEnd), cimag(tbEnd));
		}

		if (numIm > fps * lengthAnim) return(1);
	}
	return 1;
}
