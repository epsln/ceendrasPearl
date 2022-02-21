#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <complex.h>
#include <pthread.h>

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
#define WIDTH  1 * 1080 * ANTIALPOW 
#define HEIGHT 1 * 1080 * ANTIALPOW
#define BOUNDS 1 
#define RANDBOUNDS 0 + 1 * I 
#define EPSI  0.001 
#define MAXWORD 14 
#define LINE 1 
#define BITWISE 0
#define DEBUG 0


int main(){
	time_t pt;
	srand((unsigned) time(&pt));

	int numIm = 0;

	double complex ta = 0.;
	double complex tb = 0.;
	double complex tab = 0.;

	int fps = 30;
	int duration = 5;
	int lengthAnim = 10;

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

	double complex p1 = 0;
	double complex p2 = 0;


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

	p1 = randomComplex(-2 - 1. * I, 2 + 1 * I);
	p2 = randomComplex(-2 - 1. * I, 2 + 1 * I);


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

	int denum = 10;//The maximum denominator we should attain in the farray sequence
	float t = 0;

	ratio fareySeq[10 * 10];//Allocating an array for the farray sequence using the limit of its length  
	ratio fract;
	int wordLength = fract.p + fract.q;

	
	char* speWord;//Array containing the specialWord
	double complex*  fixRep; //2D array containing the fix point of all cyclic perm and inverse of the speWord
	int numFP[4] = {0}; //Number of fixed point per gens

	//makeFareySeq(denum, fareySeq);
	//makeFiboSeq(10, fareySeq);
	//makePiSeq(denum * denum, fareySeq);
	//for(int i = 1; i <= 100; i++){
	//	fareySeq[i - 1] = (ratio){2, 2*i + 1};
	//}
	//makeContinuedFraction(30, 0.1415926, fareySeq);
	dfsArgs * args = malloc(sizeof(*args));
        pthread_t threadArray[4];

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

		t   = easeInOutQuad((float)(numIm%(fps * duration) - 1), 0, 1, (float)fps * duration); 
		ta  = bezier(taBeg, p1, p2,  taEnd, t);
		tb  = bezier(tbBeg, p1, p2,  tbEnd, t);
		tab = bezier(tabBeg, p1, p2, tabEnd, t);
		tb  = easeInOutQuad((float)(numIm%(fps*duration)), tbBeg, -copysign(creal(tbBeg- tbEnd), creal(tbBeg- tbEnd)) + I*-copysign(cimag(tbBeg- tbEnd), cimag(tbBeg- tbEnd)), (float)fps * duration);
		tab = easeInOutQuad((float)(numIm%(fps*duration)), tabBeg, -copysign(creal(tabBeg- tabEnd), creal(tabBeg- tabEnd)) + I*-copysign(cimag(tabBeg- tabEnd), cimag(tabBeg- tabEnd)), (float)fps * duration);

		if (DEBUG == 1){
			printf("tab: %lf + I %lf\n", creal(tab), cimag(tab));
			printf("p/q: %lld/%lld\n", fareySeq[numIm].p, fareySeq[numIm].q);
		}

		//Compute the associated mu value...
		//fract = fareySeq[numIm + 1];
		fract = (ratio) {0, 1}; 
		//printf("p/q: %lld/%lld\n", fract.p, fract.q);
		wordLength = fract.p + fract.q;
		//getTraceFromFract(pMu, fract);
		//newtonSolver(pMu, fract); 
		//grandmaRecipe(-I*mu, 2, gens);
		grandmaRecipe(ta, 2, gens);

		//Care with the calloc !
		speWord = calloc(fract.p + fract.q, sizeof(char));
		fixRep  = calloc(4 * (fract.p + fract.q + 4), sizeof(double complex));//4 gens * (p + q number of perm + special Word abAB) 
		for (int i = 0; i < 4; i++)
			numFP[i] = 0;

			
		getSpecialWordFromFract(fract, speWord);
		computeRepetendsv2(gens, fixRep, numFP, speWord, wordLength);
		
		//printf("mu: %lf + %lf\n", creal(mu), cimag(mu));

		//Compute some generators using a recipe...
		tb = 2; 
		double complex p = -ta * tb;
		double complex q = cpow(ta, 2) + cpow(tb, 2) - 2;
		tab = (-p+csqrt(cpow(p, 2) - 4 * q))/2; 
		
		//printf("%lld/%lld\n", fareySeq[numIm].p, fareySeq[numIm].q);
		//grandmaSpecialRecipe(ta, tb, tab, gens);
		dfsArgs *args;

		//Explore depth first combination of generators...
		for (int i = 0; i < 4; i++){
			//Put me in a separate function pls
			args = malloc(sizeof(dfsArgs));
			args->gens = gens;
			args->img = pImg;
			args->numIm = numIm;
			//Need to implement a way to not hardcode the size of fixRep...
			args->specialWord = calloc(fract.p + fract.q, sizeof(char));
			args->fixRep = calloc(4 * (fract.p + fract.q + 4), sizeof(double complex));
			strcpy(args->specialWord, speWord);
			for (int i = 0; i < fract.p + fract.q; i++)
				args->specialWord[i] = speWord[i];
			
			for (int i = 0; i < 4; i++){
				args->numFP[i] = numFP[i];
			}

			for (int j = 0; j < 4; j++){
				for (int h = 0; h < numFP[j]; h++){
					args->fixRep[j * (fract.p + fract.q + 4) + h] = fixRep[j * (fract.p + fract.q + 4) + h];
				}
			}
			args->wordLength = wordLength;
			args->numBranch = modulo(4 - i, 4);
			pthread_create(&threadArray[i], NULL, computeDepthFirst, args);
		}

		//Save as bmp
		for(int i = 0; i < 4 ; i++){
			pthread_join(threadArray[i], NULL);
		}
		saveArrayAsBMP(pImg);

		//Update progress bar
		pBarAnim(numIm, fps * lengthAnim, timeArray); 

		numIm ++;
		if (numIm % (fps * duration) == 0 ){//Change target traces once we have arrived 
			taBeg = taEnd;
			tbBeg = tbEnd;
			tabBeg = tabEnd;

			taEnd  = randomComplexFixDist(taBeg, 0.25);
			tbEnd  = randomComplexFixDist(tbBeg, 0.25);
			tabEnd = randomComplexFixDist(tabBeg, 0.25);
			p1 = randomComplexFixDist(p1, 0.25);
			p2 = randomComplexFixDist(p2, 0.25);

			if (numIm >= fps * lengthAnim - fps * duration ){//loop by ending up at the begining traces
				taEnd = taInit;
				tbEnd = tbInit;
				tabEnd = tabInit;
			}
			printf("Switch ! \n");
			printf("ta : %lf %lf \n", creal(taEnd), cimag(taEnd));
			printf("tb : %lf %lf \n", creal(tbEnd), cimag(tbEnd));
			printf("Switch ! \n");
		}

		if (numIm >= fps * lengthAnim){ pthread_exit(NULL);  return(1);}//Get out of here when we're done !
		//if (fract.p == 0 && fract.q == 0){ pthread_exit(NULL);  return(1);}//Get out of here when we're done !
		//Else, we go again !
	}
	return 0;
}
