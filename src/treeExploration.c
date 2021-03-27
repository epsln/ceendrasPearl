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
#include "include/recipes.h"

void goForward(int *lev, int* tag, int* state, int FSA[19][4], double complex* word, double complex* gens){
	double complex buffWord[2][2];
	double complex buffGen[2][2];
	double complex buffOut[2][2];

	int idxGen = 0, i = 1;
	
	//printf("goForward: %d\n", *lev);
	
	for (int i = 0; i < 2; i++){
		for (int j = 0; j < 2; j++){
			buffWord[i][j] = 0;
			buffGen[i][j] = 0;
			buffOut[i][j] = 0;
		}
	}

	//printf("state[%d] = %d\n", *lev, state[*lev]);
	//printf("FSA[%d][%d] = %d\n\n", state[*lev], modulo(tag[*lev - 1] + 1, 4),  FSA[state[*lev]][modulo(tag[*lev - 1] + 1, 4)]);
	//printf("FSA[%d][%d] = %d\n", state[*lev + 1], modulo(tag[*lev + 1] - 1, 4), FSA[state[*lev + 1]][modulo(tag[*lev + 1] - 1, 4)]);
	do{
		idxGen = modulo(tag[*lev] + i, 4);
		//printf("FSA[%d][%d] == %d\n",state[*lev], idxGen, FSA[state[*lev]][idxGen]);
		i--;
	}while((FSA[state[*lev]][idxGen] == 0));
	//printf("idxGen = %d\n", idxGen);


	state[*lev + 1] = FSA[state[*lev]][idxGen];
	//printf("state[%d] = %d\n", *lev + 1, state[*lev + 1]);

	*lev = *lev + 1; //Advance one lvl

	//tag[*lev] = modulo(tag[*lev - 1] + 1, 4); //Take the rightmost turn
	tag[*lev] = idxGen; //Take the rightmost turn


	matrix3dto2D(word, buffWord, *lev - 1);
	matrix3dto2D(gens, buffGen, tag[*lev]);
	matmul(buffWord, buffGen, buffOut);

	matrix2dto3D(buffOut, word, *lev);
}

void goBackwards(int *lev){
	//printf("goBack\n");
	*lev = *lev - 1;
}

int availableTurn(int *lev, int* tag, int* state, int FSA[19][4]){
//	printf("availableTurn\n");
//	printf("FSA[%d][%d] = %d\n\n", state[*lev], modulo(tag[*lev + 1] - 1, 4), FSA[state[*lev]][modulo(tag[*lev + 1] - 1, 4)]);
	if (*lev == -1)	return 0;

	if((FSA[state[*lev]][modulo(tag[*lev + 1] - 1, 4)] == 0))
		return 0;

	//if (modulo(tag[*lev + 1] - 1,  4) == modulo(tag[*lev] + 2 , 4))
	//	return 0;
	
	else
		return 1;
	
}	

void turnForward(int *lev, int *tag, int* state, int FSA[19][4], double complex* word, double complex* gens){
	double complex buffWord[2][2];
	double complex buffGen[2][2];
	double complex buffOut[2][2];
	int idxGen = 0, i = 1;
	//printf("turnForward\n");

	//tag[*lev + 1] = modulo(tag[*lev + 1] - 1,  4);
	

	do{	
		idxGen = modulo(tag[*lev + 1] - i, 4);
	//	printf("FSA[%d][%d] == %d\n",state[*lev], idxGen + 1, FSA[state[*lev]][idxGen]);
		i--;
	}while((FSA[state[*lev]][idxGen] == 0));

//	printf("idxGen = %d\n", idxGen);
	state[*lev + 1] = FSA[state[*lev]][idxGen];
//	printf("state[%d] = %d\n", *lev + 1, state[*lev + 1]);

	tag[*lev + 1] = idxGen; 
//	printf("tag[%d] = %d\n\n", *lev + 1, tag[*lev + 1]);
	//state[*lev + 1] = modulo(tag[*lev + 1] - 1,  4) ;

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




int branchTermRepetends(int lev, int* tag, double complex fixRep[4][4], double complex* word, image_t* img){
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

	//If we hit the maximum length of word (-1 because of 0 indexing, bailout !)
	//Used only in case of point mode
	//printf("lvl: %d\n", lev);
	//printf("%lf %lf\n",creal(z0), cimag(z0));
	//printf("%lf %lf\n",creal(z1), cimag(z1));
	//printf("%lf %lf\n\n",creal(z2), cimag(z2));
	if (img->line == 1 && lev == img->maxword - 1){
		return 1;
	}

	//if the word ends with a or A, use a 4th fixed point 
	if (tag[lev] % 2 == 0){
		z3 = mobiusOnPoint(buffWord, fixRep[tag[lev]][3]);
	}
	
	//Check the distance between all points, if < epsi then draw a line from z0 to z1 to z2 to (maybe) z3	
	if ((img->line == 0 && lev == img->maxword - 1) || (cabs(z0 - z1) < img->epsi && cabs(z1 - z2) < img->epsi  && cabs(z2 - z3) < img->epsi) ){
		showMatrix(buffWord, img);
		x0 = (int) map(creal(z0), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y0 = (int) map(cimag(z0), -img->bounds, img->bounds, img->h, 0);
		x1 = (int) map(creal(z1), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y1 = (int) map(cimag(z1), -img->bounds, img->bounds, img->h, 0);
		x2 = (int) map(creal(z2), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y2 = (int) map(cimag(z2), -img->bounds, img->bounds, img->h, 0);
		x3 = (int) map(creal(z3), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		y3 = (int) map(cimag(z3), -img->bounds, img->bounds, img->h, 0);

		//TODO: Use different bounds depending on generator !
		//x0 = (int) map(creal(z0), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		//y0 = (int) map(cimag(z0), -2 * img->bounds, 0 , img->h, 0);
		//x1 = (int) map(creal(z1), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		//y1 = (int) map(cimag(z1), -2 * img->bounds, 0, img->h, 0);
		//x2 = (int) map(creal(z2), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		//y2 = (int) map(cimag(z2), -2 * img->bounds, 0, img->h, 0);
		//x3 = (int) map(creal(z3), -aspectRatio * img->bounds, aspectRatio * img->bounds, 0, img->w);
		//y3 = (int) map(cimag(z3), 0, 2 * img->bounds, img->h, 0);
		//
		//printf("%lf %lf\n",creal(z0), cimag(z0));
		//printf("%lf %lf\n",creal(z1), cimag(z1));
		//printf("%lf %lf\n",creal(z2), cimag(z2));
		line(x0, y0, x1, y1, img);	
		line(x1, y1, x2, y2, img);	
		point(x0, y0, img);
		point(x1, y1, img);
		point(x2, y2, img);

		//if the word ends with a or A, use a 4th fixed point 
		if (tag[lev] % 2 == 0){
			line(x2, y2, x3, y3, img);	
			point(x3, y3, img);
		}
		return 1;
	}
	return 0;
}


void computeDepthFirst(double complex* gens, image_t* img, int numIm){
	//TODO: clean me up
	int lev = 0;

	double complex endpt[4];//Deprecated !
	double complex begpt[4];//Deprecated !

	double complex fixRep[4][4];//Contains the fixed points of a few special combinations of gens
	
	//Array of size [maxword, 2, 2]
	//Contains all of the word computed until the current level
	double complex *word = (double complex *)calloc(img->maxword * 2 * 2, sizeof(double complex));
	//The tag array contains the index of the used gen up until current level
	int *tag  = (int *)calloc(img->maxword, sizeof(int));
	//The state array contains the state of the automaton (see pp. 360)
	//This automaton prevents us from making forbiden association that ends up in identity (like a and A) 
	int *state = (int *)calloc(img->maxword, sizeof(int));

	memset(state, 0, img->maxword * sizeof(int));
	state[0] = 1;//We start at a
	tag[0] = 0;

	//State automaton array:
	//Given the current state and the next generator, gives the next state
	//If possible, try to find an algo that generates this sort of stuff
	
	int FSA[19][4] = {
		{ 1,  2,  3,  4},//Identity 0 
		{ 1,  5,  0,  4},//a        1
		{ 7,  2,  8,  0},//b        2
		{ 0,  2,  3,  6},//A        3
		{ 9,  0, 10,  4},//B        4
		{ 7,  2, 11,  0},//ab       5
		{12,  0, 10,  4},//AB       6
		{ 1,  5,  0, 13},//ba       7
		{ 0,  2,  3, 14},//bA       8
		{ 1, 15,  0,  4},//Ba       9 
		{ 0, 16,  3,  4},//BA       10
		{ 0,  2,  3, 17},//abA      11
		{ 1, 18,  0,  4},//ABa      12
		{ 9,  0,  0,  4},//baB      13
		{ 0,  0, 10,  4},//bAB      14
		{ 7,  2,  0,  0},//Bab      15
		{ 0,  2,  8,  0},//BAb      16
		{ 0,  0, 10,  4},//abAB     17
		{ 7,  2,  0,  0} //ABab     18
	};
		
	//int FSA[5][4] = {
	//	{ 1,  2,  3,  4},//Identity 0 
	//	{ 1,  2,  0,  4},//a        1
	//	{ 1,  2,  3,  0},//b        2
	//	{ 0,  2,  3,  4},//A        3
	//	{ 1,  0,  3,  4}//B        4
	//};
	int *plev;
	plev = &lev;

	computeCycles(begpt, endpt, gens);
	computeRepetendsv2(gens, fixRep);

	matrix3dto3D(gens, word, 0, 0);

	while (!(lev == -1 && tag[0] == 1)){//See pp.148 for algo
		while(branchTermRepetends(lev, tag, fixRep, word, img) == 0){
			goForward(plev, tag, state, FSA, word, gens);	
			printWord(lev, tag, img);
		}
		do{
			goBackwards(plev);
			printWord(lev, tag, img);
		}while(!(availableTurn(plev, tag, state, FSA) == 1 || lev == -1));

		//Might need to put another exit condition check right here...
		if (lev == -1 && tag[0] == 1) break;	

		turnForward(plev, tag, state, FSA, word, gens);
		printWord(lev, tag, img);
	}
}


