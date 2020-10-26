#include <complex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "include/arraysOps.h"
#include "include/complexMath.h"
#include "include/debugTools.h"
#include "include/plot.h"

void goForward(int *lev, int* tag, double complex word[1000][2][2], double complex gens[4][2][2]){
	double complex buffWord[2][2];
	double complex buffGen[2][2];
	double complex buffOut[2][2];
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			buffWord[i][j] = 0;
			buffGen[i][j] = 0;
			buffOut[i][j] = 0;
		}
	}

	*lev = *lev + 1; //Advance one lvl

	tag[*lev] = modulo(tag[*lev - 1] + 1, 4); //Take the rightmost turn


	matrix3dto2D(word, buffWord, *lev - 1);
	matrix3dto2D(gens, buffGen, tag[*lev]);
	matmul(buffWord, buffGen, buffOut);

	matrix2dto3D(buffOut, word, *lev);
}

void goBackwards(int *lev){
	*lev = *lev - 1;
}

int availableTurn(int *lev, int* tag){

	if (modulo(tag[*lev + 1] - 1,  4) == modulo(tag[*lev] + 2 , 4)){
		return 0;
	}
	else{
		return 1;
	}
}	

void turnForward(int *lev, int tag[100000], double complex word[100000][2][2], double complex gens[4][2][2]){
	double complex buffWord[2][2];
	double complex buffGen[2][2];
	double complex buffOut[2][2];
	tag[*lev + 1] = modulo(tag[*lev + 1] - 1,  4);
	if (*lev == -1)
		matrix3dto3D(gens, word, 0, tag[0]);
	else{
		matrix3dto2D(word, buffWord, *lev);
		matrix3dto2D(gens, buffGen, tag[*lev + 1]);
		matmul(buffWord, buffGen, buffOut);
		matrix2dto3D(buffOut, word, *lev + 1);
	}
	*lev = *lev + 1;
}



int branchTermEpsi(double* PARAMS, double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex word[1000][2][2], float*** imgArr){
	//Basic branch term using only the distance between an older point and a new point
	//See pp. 185
	//TODO: Find better way to pass parameters between functions
	int LEVMAX    = (int) PARAMS[0];
	double EPSI   = PARAMS[1];
	double BOUNDS = PARAMS[2];
	int WIDTH     = (int) PARAMS[3];
	int HEIGHT    = (int) PARAMS[4];
	int LINE      = (int) PARAMS[5];

	double complex buffWord[2][2];
	matrix3dto2D(word, buffWord, lev);
	double complex newPoint = mobiusOnPoint(buffWord, endpt[tag[lev]]);
	showMatrix(buffWord, PARAMS);
	if (lev == LEVMAX || cabs(newPoint - *oldPoint) < EPSI){
		int x0 = (int) map(creal(newPoint), -BOUNDS, BOUNDS, 0, WIDTH);
		int y0 = (int) map(cimag(newPoint), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x1 = (int) map(creal(*oldPoint), -BOUNDS, BOUNDS, 0, WIDTH);
		int y1 = (int) map(cimag(*oldPoint), -BOUNDS, BOUNDS, HEIGHT, 0);
		if (*oldPoint != -1000){
			if (checkBoundaries(x0, y0, WIDTH, HEIGHT) == 1 && checkBoundaries(x1, y1, WIDTH, HEIGHT) == 1){
				line(x0, y0, x1, y1, imgArr, LINE, WIDTH, HEIGHT);	
				point(x0, y0, imgArr, WIDTH, HEIGHT);
				point(x1, y1, imgArr, WIDTH, HEIGHT);
			}
		}
		*oldPoint = newPoint;
		return 1;
	}
	return 0;
}

int branchTermRepetends(double* PARAMS, double complex* oldPoint, int lev, int* tag, double complex fixRep[4][4], double complex word[1000][2][2], float*** imgArr){
	//Branch termination check based on the repetends methods
	int LEVMAX    = (int) PARAMS[0];
	double EPSI   = PARAMS[1];
	double BOUNDS = PARAMS[2];
	int WIDTH     = (int) PARAMS[3];
	int HEIGHT    = (int) PARAMS[4];
	int LINE      = (int) PARAMS[5];
	int DEBUG     = (int) PARAMS[5];
	double complex buffWord[2][2];
	matrix3dto2D(word, buffWord, lev);

	double complex z0 = mobiusOnPoint(buffWord, fixRep[tag[lev]][0]);
	double complex z1 = mobiusOnPoint(buffWord, fixRep[tag[lev]][1]);
	double complex z2 = mobiusOnPoint(buffWord, fixRep[tag[lev]][2]);
	double complex z3 = z2;

	if (tag[lev] % 2 == 0){
		z3 = mobiusOnPoint(buffWord, fixRep[tag[lev]][3]);
	}

	if (lev == LEVMAX || (cabs(z0 - z1) < EPSI && cabs(z1 - z2) < EPSI  && cabs(z2 - z3) )){
		//showMatrix(buffWord);
		int x0 = (int) map(creal(z0), -BOUNDS, BOUNDS, 0, WIDTH);
		int y0 = (int) map(cimag(z0), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x1 = (int) map(creal(z1), -BOUNDS, BOUNDS, 0, WIDTH);
		int y1 = (int) map(cimag(z1), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x2 = (int) map(creal(z2), -BOUNDS, BOUNDS, 0, WIDTH);
		int y2 = (int) map(cimag(z2), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x3 = (int) map(creal(z2), -BOUNDS, BOUNDS, 0, WIDTH);
		int y3 = (int) map(cimag(z2), -BOUNDS, BOUNDS, HEIGHT, 0);
		line(x0, y0, x1, y1, imgArr, LINE, WIDTH, HEIGHT);	
		line(x1, y1, x2, y2, imgArr, LINE, WIDTH, HEIGHT);	
		point(x0, y0, imgArr, WIDTH, HEIGHT);
		point(x1, y1, imgArr, WIDTH, HEIGHT);
		if (tag[lev] % 2 == 0){
			line(x2, y2, x3, y3, imgArr, LINE, WIDTH, HEIGHT);	
			point(x0, y0, imgArr, WIDTH, HEIGHT);
		}
		*oldPoint = z2;
		return 1;
	}
	return 0;
}


void computeDepthFirst(double* PARAMS, double complex ta, double complex tb, double complex tab, float*** imgArr, int numIm){
	int lev = 0;
	char filename[100] = "out/img_";
	char imageNum[5];  
	sprintf(imageNum, "%d", numIm);
	strcat(filename, imageNum);
	strcat(filename, ".bmp");
	printf("%s\n", filename);


	double complex oldPoint = 0.;

	double complex endpt[4];
	double complex begpt[4];
	double complex gens[4][2][2];
	double complex fixRep[4][4];
	double complex word[1000][2][2];
	int tag[1000];
	int *plev;
	plev = &lev;
	double complex *poldP = &oldPoint;

	for (int i = 0; i < 1000; i++){
		tag[i] = 0;
	}

//	grandmaRecipe(ta, tb, gens);
	grandmaSpecialRecipe(ta, tb, tab, gens);
	printf("a = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[0][0][0]),cimag(gens[0][0][0]), creal(gens[0][1][0]), cimag(gens[0][1][0]), creal(gens[0][0][1]),cimag(gens[0][0][1]), creal(gens[0][1][1]),cimag(gens[0][1][1]));
	printf("b = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[1][0][0]),cimag(gens[1][0][0]), creal(gens[1][1][0]), cimag(gens[1][1][0]), creal(gens[1][0][1]),cimag(gens[1][0][1]), creal(gens[1][1][1]),cimag(gens[1][1][1]));
	printf("A = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[2][0][0]),cimag(gens[2][0][0]), creal(gens[2][1][0]), cimag(gens[2][1][0]), creal(gens[2][0][1]),cimag(gens[2][0][1]), creal(gens[2][1][1]),cimag(gens[2][1][1]));
	printf("B = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[3][0][0]),cimag(gens[3][0][0]), creal(gens[3][1][0]), cimag(gens[3][1][0]), creal(gens[3][0][1]),cimag(gens[3][0][1]), creal(gens[3][1][1]),cimag(gens[3][1][1]));

	computeCycles(begpt, endpt, gens);
	computeRepetendsv2(gens, fixRep);

	matrix3dto3D(gens, word, 0, 0);
	printf("word[0] = [[%lf + i %lf, %lf + i%lf],\n[%lf + i %lf, %lf + i %lf ]]\n", creal(word[0][0][0]),cimag(word[0][0][0]), creal(word[0][1][0]), cimag(word[0][1][0]), creal(word[0][0][1]),cimag(word[0][0][1]), creal(word[0][1][1]),cimag(word[0][1][1]));

	oldPoint = begpt[3];
	while (!(lev == -1 && tag[0] == 1)){//See pp.148 for algo
		while(branchTermRepetends(PARAMS, poldP, lev, tag, fixRep, word, imgArr) == 0){
			goForward(plev, tag, word, gens);	
			printWord(lev, tag, PARAMS);
		}
		do{
			goBackwards(plev);
			printWord(lev, tag, PARAMS);
		}while(!(availableTurn(plev, tag) == 1 || lev == -1));

		//Might need to put another exit condition check right here...
		if (lev == -1 && tag[0] == 1) break;	
		turnForward(plev, tag, word, gens);
		printWord(lev, tag, PARAMS);
	}
	saveArrayAsBMP(imgArr,filename, PARAMS);

}


