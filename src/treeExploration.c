#include <complex.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "include/arraysOps.h"
#include "include/complexMath.h"
#include "include/debugTools.h"
#include "include/plot.h"
#include "include/treeExploration.h"

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
	if (*lev == -1)	return 0;

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



int branchTermEpsi(double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex word[1000][2][2], image_t* img){
	//Basic branch term using only the distance between an older point and a new point
	//See pp. 185

	float aspectRatio = img->w/(float)img->h;

	double complex buffWord[2][2];
	matrix3dto2D(word, buffWord, lev);

	double complex newPoint = mobiusOnPoint(buffWord, endpt[tag[lev]]);

	showMatrix(buffWord, img);
	int x0, x1;
	int y0, y1;

	if (lev == img->levmax || cabs(newPoint - *oldPoint) < img->epsi){
		x0 = (int) map(creal(newPoint), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y0 = (int) map(cimag(newPoint), -img->bounds, img->bounds, img->h, 0);
		x1 = (int) map(creal(*oldPoint), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y1 = (int) map(cimag(*oldPoint), -img->bounds, img->bounds, img->h, 0);

		line(x0, y0, x1, y1, img);	
		point(x0, y0, img);
		point(x1, y1, img);
		*oldPoint = newPoint;

		return 1;
	}

	return 0;
}

int branchTermRepetends(double complex* oldPoint, int lev, int* tag, double complex fixRep[4][4], double complex word[1000][2][2], image_t* img){
	//Branch termination check based on the repetends methods
	float aspectRatio = img->w/(float)img->h;
	double complex buffWord[2][2];
	matrix3dto2D(word, buffWord, lev);

	double complex z0 = mobiusOnPoint(buffWord, fixRep[tag[lev]][0]);
	double complex z1 = mobiusOnPoint(buffWord, fixRep[tag[lev]][1]);
	double complex z2 = mobiusOnPoint(buffWord, fixRep[tag[lev]][2]);
	double complex z3 = z2;

	int x0, x1, x2, x3;
	int y0, y1, y2, y3;

	if (tag[lev] % 2 == 0){
		z3 = mobiusOnPoint(buffWord, fixRep[tag[lev]][3]);
	}

	if (lev == img->levmax || (cabs(z0 - z1) < img->epsi && cabs(z1 - z2) < img->epsi  && cabs(z2 - z3) )){
		showMatrix(buffWord, img);
		x0 = (int) map(creal(z0), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y0 = (int) map(cimag(z0), -img->bounds, img->bounds, img->h, 0);
		x1 = (int) map(creal(z1), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y1 = (int) map(cimag(z1), -img->bounds, img->bounds, img->h, 0);
		x2 = (int) map(creal(z2), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y2 = (int) map(cimag(z2), -img->bounds, img->bounds, img->h, 0);
		x3 = (int) map(creal(z3), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y3 = (int) map(cimag(z3), -img->bounds, img->bounds, img->h, 0);
		line(x0, y0, x1, y1, img);	
		line(x1, y1, x2, y2, img);	
		point(x0, y0, img);
		point(x1, y1, img);
		point(x2, y2, img);
		if (tag[lev] % 2 == 0){
			line(x2, y2, x3, y3, img);	
			point(x3, y3, img);
		}
		*oldPoint = z2;
		return 1;
	}
	return 0;
}


void computeDepthFirst(double complex ta, double complex tb, double complex tab, image_t* img, int numIm){
	int lev = 0;
	
	double complex oldPoint = 0.;

	double complex endpt[4];
	double complex begpt[4];
	double complex gens[4][2][2];
	double complex fixRep[4][4];
	double complex word[1000][2][2];
	int tag[1000] = {0};
	int *plev;
	plev = &lev;
	double complex *poldP = &oldPoint;

	grandmaRecipe(ta, tb, gens);

//	grandmaSpecialRecipe(ta, tb, tab, gens);
//	printf("a = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[0][0][0]),cimag(gens[0][0][0]), creal(gens[0][1][0]), cimag(gens[0][1][0]), creal(gens[0][0][1]),cimag(gens[0][0][1]), creal(gens[0][1][1]),cimag(gens[0][1][1]));
//	printf("b = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[1][0][0]),cimag(gens[1][0][0]), creal(gens[1][1][0]), cimag(gens[1][1][0]), creal(gens[1][0][1]),cimag(gens[1][0][1]), creal(gens[1][1][1]),cimag(gens[1][1][1]));
//	printf("A = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[2][0][0]),cimag(gens[2][0][0]), creal(gens[2][1][0]), cimag(gens[2][1][0]), creal(gens[2][0][1]),cimag(gens[2][0][1]), creal(gens[2][1][1]),cimag(gens[2][1][1]));
//	printf("B = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[3][0][0]),cimag(gens[3][0][0]), creal(gens[3][1][0]), cimag(gens[3][1][0]), creal(gens[3][0][1]),cimag(gens[3][0][1]), creal(gens[3][1][1]),cimag(gens[3][1][1]));

	computeCycles(begpt, endpt, gens);
	computeRepetendsv2(gens, fixRep);

	matrix3dto3D(gens, word, 0, 0);
	

	oldPoint = begpt[3];
	while (!(lev == -1 && tag[0] == 1)){//See pp.148 for algo
		while(branchTermRepetends(poldP, lev, tag, fixRep, word, img) == 0){
			goForward(plev, tag, word, gens);	
			printWord(lev, tag, img);
		}
		do{
			goBackwards(plev);
			printWord(lev, tag, img);
		}while(!(availableTurn(plev, tag) == 1 || lev == -1));

		//Might need to put another exit condition check right here...
		if (lev == -1 && tag[0] == 1) break;	
		turnForward(plev, tag, word, gens);
		printWord(lev, tag, img);
	}
}


