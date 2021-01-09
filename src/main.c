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

#define ANTIALPOW 4
#define WIDTH  1000 * ANTIALPOW 
#define HEIGHT 1000 * ANTIALPOW
#define BOUNDS 1 
#define RANDBOUNDS 0 + 1 * I 
#define EPSI  0.01 
#define MAXWORD 1000 
#define LINE 1 
#define BITWISE 0
#define DEBUG 0


int main(){
	time_t pt;
	srand((unsigned) time(&pt));

	int fps = 10;
	int duration = 3;
	int lengthAnim = 1;

	double complex ta = 0.;
	double complex tb = 0.;
	double complex tab = 0.;

	double complex taBeg = 0.;
	double complex tbBeg = 0.;
	double complex tabBeg = 0.;

	double complex taEnd =  0.;
	double complex tbEnd = 0.;
	double complex tabEnd =  0.;

	double complex taInit = 0.;
	double complex tbInit = 0.;
	double complex tabInit = 0.;

	int numIm = 0;


	taBeg  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
	tbBeg  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
	tabBeg = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);

	taEnd  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
	tbEnd  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
	tabEnd = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);

	tbBeg  = 2;
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
	makeFiboSeq(1000, fareySeq);

	while(1){
		srand((unsigned) time(&pt));
		sprintf(imageNum, "%d", numIm);
		strcat(prefix, imageNum);
		strcat(prefix, ".bmp\0");
		strcpy(pImg->filename, prefix);
		strcpy(prefix, "out/img_");


		//Create a filename for the image based on the number of image processed
		//TODO: Move this to its own function :)

		//printf("Image: %s\n\n", pImg->filename);

		//Here, we interpolate between two traces using an easing function
		//ta = InOutQuadComplex((float)(numIm%(fps*duration)), taBeg, -copysign(creal(taBeg- taEnd), creal(taBeg- taEnd)) + I*-copysign(cimag(taBeg- taEnd), cimag(taBeg- taEnd)), (float)fps * duration); tb = InOutQuadComplex((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		//tb = InOutQuadComplex((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		//tab = InOutQuadComplex((float)(numIm%(fps*duration)), tabBeg, -copysign(creal(tabBeg- tabEnd), creal(tabBeg- tabEnd)) + I*-copysign(cimag(tabBeg- tabEnd), cimag(tabBeg- tabEnd)), (float)fps * duration);

		if (DEBUG == 1){
			printf("ta:  %lf + I %lf\n", creal(ta), cimag(ta));
			printf("tb:  %lf + I %lf\n", creal(tb), cimag(tb));
			printf("tab: %lf + I %lf\n", creal(tab), cimag(tab));
		}

		printf("p/q: %lld/%lld\n", fareySeq[numIm].p, fareySeq[numIm].q);
		//Compute the associated mu value...
		newtonSolver(pMu, fareySeq[numIm + 2]);
		printf("mu: %lf + %lf\n", creal(mu), cimag(mu));

		//Compute some generators using a recipe...
		//maskitRecipe(mu, gens);
		//grandmaRecipe(2, 2, gens);
		grandmaRecipe(-I*mu, 2, gens);

		//Explore depth first combination of generators...
		computeDepthFirst(gens, pImg, numIm);
		saveArrayAsBMP(pImg);
		numIm ++;
		//if (numIm % (fps * duration) == 0 ){//Change target traces once we have arrived 
		//	taBeg = taEnd;
		//	tbBeg = tbEnd;
		//	tabBeg = tabEnd;

		//	taEnd  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
		//	tbEnd  = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);
		//	tabEnd = randomComplex(-3 - 1.5 * I, 3 + 1.5 * I);


		//	if (numIm >= fps * lengthAnim - fps * duration ){//loop by ending up at the begining traces
		//		taEnd = taInit;
		//		tbEnd = tbInit;
		//		tabEnd = tabInit;
		//	}
		//}

		if (numIm > 69) return(1);//Get out of here when we're done !
		//if (fareySeq[numIm].p == 0 && fareySeq[numIm].q == 0) return(1);//Get out of here when we're done !
		//if (fareySeq[numIm].p == 0 && fareySeq[numIm].q == 0) return(1);//Get out of here when we're done !
		//if (numIm >= fps * lengthAnim) return(1);//Get out of here when we're done !
		//Else, we go again !
	}
	return 0;
}
