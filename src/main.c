#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <complex.h>

#include "include/complexMath.h"
#include "include/plot.h"
#include "include/arraysOps.h"
#include "include/treeExploration.h"
#include "include/debugTools.h"
#include "include/easing.h"
#include "include/recipes.h"
#include "include/accidents.h"
#include "include/progressBar.h"

#define ANTIALPOW 4
#define WIDTH  2000 * ANTIALPOW 
#define HEIGHT 2000 * ANTIALPOW
#define BOUNDS 1 
#define RANDBOUNDS 0 + 1 * I 
#define EPSI  0.00001 
#define MAXWORD 100 
#define LINE 1 
#define BITWISE 1
#define DEBUG 0


int main(){
	time_t pt;
	srand((unsigned) time(&pt));

	int numIm = 0;

	double complex ta = 0.;
	double complex tb = 0.;
	double complex tab = 0.;

	int fps = 60;
	int duration = 5;
	int lengthAnim = 30;

	//Using number of clock ticks to estimate time
	//Not using a time_t timestamp in case of subsecond compute time
	//Use a rolling average of the past 10 times to get something semi accurate
	double timeArray[11];
	


	double complex taBeg = 0.;
	double complex tbBeg = 0.;
	double complex tabBeg = 0.;

	double complex taEnd =  0.;
	double complex tbEnd = 0.;
	double complex tabEnd =  0.;

	double complex taInit = 0.;
	double complex tbInit = 0.;
	double complex tabInit = 0.;


	taBeg  = randomComplex(-2 - 1. * I, 2 + 1 * I);
	tbBeg  = randomComplex(-2 - 1. * I, 2 + 1 * I);
	tabBeg = randomComplex(-2 - 1. * I, 2 + 1 * I);

	taEnd  = randomComplex(-2 - 1. * I, 2 + 1 * I);
	tbEnd  = randomComplex(-2 - 1. * I, 2 + 1 * I);
	tabEnd = randomComplex(-2 - 1. * I, 2 + 1 * I);
	tbEnd  = 2; 

	taInit = taBeg;
	tbInit = tbBeg;
	tabInit = tabBeg;


	double complex* gens = (double complex*)calloc(4*2*2, sizeof(double complex));

	//TODO: Move this portion to its own file :)
	image_t img;
	image_t* pImg = &img;

	pImg->w       = WIDTH;
	pImg->h       = HEIGHT;
	pImg->bounds  = BOUNDS;
	pImg->epsi    = EPSI;
	pImg->line    = LINE;
	pImg->maxword = MAXWORD;
	pImg->antialiasingPow = ANTIALPOW;
	pImg->debug  = DEBUG;
	pImg->bitwise= BITWISE;
	pImg->filename = malloc(256* sizeof(char));

	pImg->pointArr = NULL;
	pImg->bitArray = NULL;
	pImg->pointArr = (int*)calloc(pImg->w*pImg->h, sizeof(int));
	unsigned int allocSize = (ceil(pImg->w/64.0))* pImg->h;
	pImg->bitArray = (unsigned long long int*)calloc(allocSize, sizeof(unsigned long long int));

	if (pImg->pointArr == NULL){
		printf("Could not allocate memory for the image array !\nExiting...\n");
		exit(-1);
	}

	char prefix[100] = "out/img_";
	char imageNum[6];  

	double complex mu = 2*I;
	double complex *pMu = &mu; 

	int denum = 25;//The maximum denominator we should attain in the farray sequence

	//ratio *fareySeq = (ratio * ) malloc(denum*denum);//Allocating an array for the farray sequence using the limit of its length  
	ratio fareySeq[denum*denum];//Allocating an array for the farray sequence using the limit of its length  

	//makeFareySeq(denum, fareySeq);
	//makeFiboSeq(1000, fareySeq);
	//makePiSeq(denum * denum, fareySeq);
	//for(int i = 1; i <= 30; i++){
	//	fareySeq[i - 1] = (ratio){1, i};
	//}
	//makeContinuedFraction(10, (sqrt(35) - 5)/10.0, fareySeq);

	while(1){
		//Create a filename for the image based on the number of image processed
		//TODO: Move this to its own function :)
		srand((unsigned) time(&pt));
		sprintf(imageNum, "%d", numIm);
		strcat(prefix, imageNum);
		strcat(prefix, ".bmp\0");
		strcpy(pImg->filename, prefix);
		strcpy(prefix, "out/img_");

		//printf("Image: %s\n\n", pImg->filename);

		//Here, we interpolate between two traces using an easing function
		ta = InOutQuadComplex((float)(numIm%(fps*duration)), taBeg, -copysign(creal(taBeg- taEnd), creal(taBeg- taEnd)) + I*-copysign(cimag(taBeg- taEnd), cimag(taBeg- taEnd)), (float)fps * duration); tb = InOutQuadComplex((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		tb = InOutQuadComplex((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		tab = InOutQuadComplex((float)(numIm%(fps*duration)), tabBeg, -copysign(creal(tabBeg- tabEnd), creal(tabBeg- tabEnd)) + I*-copysign(cimag(tabBeg- tabEnd), cimag(tabBeg- tabEnd)), (float)fps * duration);

		if (DEBUG == 1){
			printf("tab: %lf + I %lf\n", creal(tab), cimag(tab));
			printf("p/q: %lld/%lld\n", fareySeq[numIm].p, fareySeq[numIm].q);
		}
		printf("ta:  %lf + I %lf\n", creal(ta), cimag(ta));
		printf("tab:  %lf + I %lf\n", creal(tab), cimag(tab));

		//Compute the associated mu value...
		//newtonSolver(pMu, fareySeq[numIm]);
		//printf("mu: %lf + %lf\n", creal(mu), cimag(mu));

		//Compute some generators using a recipe...
		//grandmaRecipe(-I*mu, 3, gens);
		grandmaRecipe(2, 2, gens);
		//grandmaSpecialRecipe(2, ta, tab, gens);

		//Explore depth first combination of generators...
		computeDepthFirst(gens, pImg, numIm);

		//Save as bmp
		saveArrayAsBMP(pImg);

		//Update progress bar
		pBarAnim(numIm, fps * lengthAnim, timeArray); 
		exit(-1);

		numIm ++;
		if (numIm % (fps * duration) == 0 ){//Change target traces once we have arrived 
			taBeg = taEnd;
			tbBeg = tbEnd;
			tabBeg = tabEnd;

			taEnd  = randomComplex(-2 - 1 * I, 2 + 1 * I);
			tbEnd  = randomComplex(-2 - 1 * I, 2 + 1 * I);
			tabEnd = randomComplex(-2 - 1 * I, 2 + 1 * I);

			tbEnd = 2;

			if (numIm >= fps * lengthAnim - fps * duration ){//loop by ending up at the begining traces
				taEnd = taInit;
				tbEnd = tbInit;
				tabEnd = tabInit;
			}
		}

		//if (numIm > 30) return(1);//Get out of here when we're done !
		//if (fareySeq[numIm].p == 0 && fareySeq[numIm].q == 0) return(1);//Get out of here when we're done !
		//if (fareySeq[numIm].p == 0 && fareySeq[numIm].q == 0) return(1);//Get out of here when we're done !
		if (numIm >= fps * lengthAnim) return(1);//Get out of here when we're done !
		//Else, we go again !
	}
	return 0;
}
