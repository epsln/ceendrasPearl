#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include <complex.h>

#include "include/plot.h"
#include "include/complexMath.h"

#define SIZEARR 10000
#define HEIGHT 1000 
#define WIDTH  1000 
#define BOUNDS 2 
#define EPSI  0.001 
#define LINE 0 

int LEVMAX = 8;  

#define DEBUG 0

struct circle{
	complex c;
	double r;
} circle;

double map(double n,double  start1,double  stop1,double  start2,double  stop2){//Classic mapping formula
	return ((n-start1)/(stop1-start1))*(stop2-start2)+start2;
}

void showMatrix(double complex mat[2][2]){
	if (DEBUG == 1)
		printf("[[%lf + i %lf, %lf + i%lf],\n[%lf + i %lf, %lf + i %lf ]]\n", creal(mat[0][0]),cimag(mat[0][0]), creal(mat[1][0]), cimag(mat[1][0]), creal(mat[0][1]),cimag(mat[0][1]), creal(mat[1][1]),cimag(mat[1][1]));

}



void matrix3dto3D(double complex matrix3d[1000][2][2], double complex destination[1000][2][2], int i, int k){
	destination[i][0][0] = matrix3d[k][0][0];
	destination[i][0][1] = matrix3d[k][0][1];
	destination[i][1][0] = matrix3d[k][1][0];
	destination[i][1][1] = matrix3d[k][1][1];
}

void matrix3dto2D(double complex matrix3d[1000][2][2], double complex destination[2][2], int i){
	destination[0][0] = matrix3d[i][0][0];
	destination[0][1] = matrix3d[i][0][1];
	destination[1][0] = matrix3d[i][1][0];
	destination[1][1] = matrix3d[i][1][1];
}

void matrix2dto3D(double complex matrix2d[2][2], double complex destination[1000][2][2], int i){
	destination[i][0][0] = matrix2d[0][0];
	destination[i][0][1] = matrix2d[0][1];
	destination[i][1][0] = matrix2d[1][0];
	destination[i][1][1] = matrix2d[1][1];
}



void printWord(int lev, int* tag){
	//	sleep(1);
	if (DEBUG == 1){
		for (int i = 0; i <= lev; i++){
			if (tag[i] == 0)
				printf("a");
			else if ( tag[i] == 1)
				printf("b");
			else if (tag[i] == 2)
				printf("A");
			else if (tag[i] == 3)
				printf("B");
			else
				printf("%d", tag[i]);
		}
		printf("\n");
	}
}

void computeRepetends(double complex gens[4][2][2], double complex fixRep[4][3]){//See pp.218
	double complex buff_gen_a[2][2];
	double complex buff_gen_b[2][2];
	double complex buff_gen_A[2][2];
	double complex buff_gen_B[2][2];
	double complex buff_out0[2][2];
	double complex buff_out1[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2D(gens, buff_gen_a, 0);
	matrix3dto2D(gens, buff_gen_b, 1);
	matrix3dto2D(gens, buff_gen_A, 2);
	matrix3dto2D(gens, buff_gen_B, 3);

	//bAba
	matmul(buff_gen_b, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][0] = fix(buff_out0);

	//a
	fixRep[0][1] = fix(buff_gen_a);

	//BAba
	matmul(buff_gen_B, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][2] = fix(buff_out0);

	//ABab
	matmul(buff_gen_A, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][0] = fix(buff_out0);

	//b
	fixRep[1][1] = fix(buff_gen_b);

	//aBAb
	matmul(buff_gen_a, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][2] = fix(buff_out0);

	//BabA
	matmul(buff_gen_B, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][0] = fix(buff_out0);

	//A
	fixRep[2][1] = fix(buff_gen_A);

	//baBA
	matmul(buff_gen_b, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][2] = fix(buff_out0);

	//abAB
	matmul(buff_gen_a, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][0] = fix(buff_out0);

	//B
	fixRep[3][1] = fix(buff_gen_B);

	//AbaB
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][2] = fix(buff_out0);
}


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
	//printf("lev %d, tag[%d] = %d\nbuffWord:\n", *lev, *lev, tag[*lev]);
	//showMatrix(buffWord);
	//showMatrix(buffGen);
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

void turnForward(int *lev, int tag[1000], double complex word[1000][2][2], double complex gens[4][2][2]){
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



int branchTerm(double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex fixRep[4][3], double complex word[1000][2][2], float*** imgArr){
	double complex buffWord[2][2];
	matrix3dto2D(word, buffWord, lev);
	double complex newPoint = mobiusOnPoint(buffWord, endpt[tag[lev]]);
	
	double complex z0 = mobiusOnPoint(buffWord, fixRep[tag[lev]][0]);
	double complex z1 = mobiusOnPoint(buffWord, fixRep[tag[lev]][1]);
	double complex z2 = mobiusOnPoint(buffWord, fixRep[tag[lev]][2]);

	//if (lev == LEVMAX || cabs(newPoint - *oldPoint) < EPSI){
	if (lev == LEVMAX || cabs(z0 - z1) < EPSI && cabs(z1 - z2) < EPSI ){
		//showMatrix(buffWord);
		//int x0 = (int) map(creal(*oldPoint), -BOUNDS, BOUNDS, 0, WIDTH);
		//int y0 = (int) map(cimag(*oldPoint), -BOUNDS, BOUNDS, HEIGHT, 0);
		//int x1 = (int) map(creal(newPoint), -BOUNDS, BOUNDS, 0, WIDTH);
		//int y1 = (int) map(cimag(newPoint), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x0 = (int) map(creal(z0), -BOUNDS, BOUNDS, 0, WIDTH);
		int y0 = (int) map(cimag(z0), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x1 = (int) map(creal(z1), -BOUNDS, BOUNDS, 0, WIDTH);
		int y1 = (int) map(cimag(z1), -BOUNDS, BOUNDS, HEIGHT, 0);
		int x2 = (int) map(creal(z2), -BOUNDS, BOUNDS, 0, WIDTH);
		int y2 = (int) map(cimag(z2), -BOUNDS, BOUNDS, HEIGHT, 0);
		if (*oldPoint != -1000){

		//	if (LINE == 1) line(x0, y0, x1, y1, imgArr);	
			line(x0, y0, x1, y1, imgArr, LINE, WIDTH, HEIGHT);	
			line(x1, y1, x2, y2, imgArr, LINE, WIDTH, HEIGHT);	
			if (checkBoundaries(x0, y0, WIDTH, HEIGHT) == 1 && checkBoundaries(x1, y1, WIDTH, HEIGHT) == 1 && checkBoundaries(x2, y2, WIDTH, HEIGHT) == 1){
			imgArr[x0][y0][0] = 255;
			imgArr[x0][y0][1] = 255;
			imgArr[x0][y0][2] = 255;
			imgArr[x1][y1][0] = 255;
			imgArr[x1][y1][1] = 255;
			imgArr[x1][y1][2] = 255;
			imgArr[x2][y2][0] = 255;
			imgArr[x2][y2][1] = 255;
			imgArr[x2][y2][2] = 255;
			}
		}
		*oldPoint = newPoint;
		return 1;
	}
	else
		return 0;
}

void computeDepthFirst(double complex ta, double complex tb, float*** imgArr, int numIm){
	int lev = 0;
	char filename[100] = "out/img_";
	char imageNum[5];  
	sprintf(imageNum, "%d", numIm);
	strcat(filename, imageNum);
	strcat(filename, ".bmp");
	printf("%s\n", filename);


	double complex oldPoint = -1000;

	double complex gens[4][2][2];
	double complex endpt[4];
	double complex fixRep[4][3];
	double complex word[1000][2][2];
	double complex group[1000][2][2];
	int tag[1000];
	int *plev;
	plev = &lev;
	double complex *poldP = &oldPoint;

	for (int i = 0; i < 1000; i++){
		tag[i] = 0;
	}
	/*
	   double complex*** word;
	   double complex*** group;
	   tag   = malloc(sizeof(int) * SIZEARR); 
	   word  = malloc(sizeof(double complex**) * SIZEARR); 
	   group = malloc(sizeof(double complex**) * SIZEARR); 
	   for (int i = 0; i < SIZEARR; i++){
	   word[i]  = (double complex**) malloc( 2 * sizeof(double complex*));
	   group[i] = (double complex**) malloc( 2 * sizeof(double complex*));
	   for (int j = 0; j < 2; j++){
	   word[i][j] = (double complex*) malloc( 2 * sizeof(double complex));
	   group[i][j] = (double complex*) malloc( 2 * sizeof(double complex));
	   }
	   }
	   */

	grandmaRecipe(ta, tb, gens);
	printf("a = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[0][0][0]),cimag(gens[0][0][0]), creal(gens[0][1][0]), cimag(gens[0][1][0]), creal(gens[0][0][1]),cimag(gens[0][0][1]), creal(gens[0][1][1]),cimag(gens[0][1][1]));
	printf("b = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[1][0][0]),cimag(gens[1][0][0]), creal(gens[1][1][0]), cimag(gens[1][1][0]), creal(gens[1][0][1]),cimag(gens[1][0][1]), creal(gens[1][1][1]),cimag(gens[1][1][1]));
	printf("A = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[2][0][0]),cimag(gens[2][0][0]), creal(gens[2][1][0]), cimag(gens[2][1][0]), creal(gens[2][0][1]),cimag(gens[2][0][1]), creal(gens[2][1][1]),cimag(gens[2][1][1]));
	printf("B = [[%lf + i %lf, %lf + i %lf],\n     [%lf + i %lf, %lf + i %lf ]]\n\n", creal(gens[3][0][0]),cimag(gens[3][0][0]), creal(gens[3][1][0]), cimag(gens[3][1][0]), creal(gens[3][0][1]),cimag(gens[3][0][1]), creal(gens[3][1][1]),cimag(gens[3][1][1]));
	double complex buff_gen0[2][2];
	double complex buff_gen1[2][2];
	double complex buff_gen2[2][2];
	double complex buff_gen3[2][2];
	double complex buff_out0[2][2];
	double complex buff_out1[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2D(gens, buff_gen0, 0);
	matrix3dto2D(gens, buff_gen1, 1);
	matrix3dto2D(gens, buff_gen2, 2);
	matrix3dto2D(gens, buff_gen3, 3);

	//Compute the fix points of all right most turns 
	matmul(buff_gen3, buff_gen2, buff_out0); 
	matmul(buff_out0, buff_gen1, buff_out1); 
	matmul(buff_out1, buff_gen0, buff_out0); 
	endpt[0] = fix(buff_out0);

	matmul(buff_gen0, buff_gen3, buff_out0); 
	matmul(buff_out0, buff_gen2, buff_out1); 
	matmul(buff_out1, buff_gen1, buff_out0); 
	endpt[1] = fix(buff_out0);

	matmul(buff_gen1, buff_gen0, buff_out0); 
	matmul(buff_out0, buff_gen3, buff_out1); 
	matmul(buff_out1, buff_gen2, buff_out0); 
	endpt[2] = fix(buff_out0);

	matmul(buff_gen2, buff_gen1, buff_out0); 
	matmul(buff_out0, buff_gen0, buff_out1); 
	matmul(buff_out1, buff_gen3, buff_out0); 
	endpt[3] = fix(buff_out0);

	for (int i = 0; i < 4; i++){
		printf("endpt[%d] = %lf + i %lf)\n",i,  creal(endpt[i]), cimag(endpt[i]));
	}

	computeRepetends(gens, fixRep);

	matrix3dto3D(gens, word, 0, 0);
	printf("word[0] = [[%lf + i %lf, %lf + i%lf],\n[%lf + i %lf, %lf + i %lf ]]\n", creal(word[0][0][0]),cimag(word[0][0][0]), creal(word[0][1][0]), cimag(word[0][1][0]), creal(word[0][0][1]),cimag(word[0][0][1]), creal(word[0][1][1]),cimag(word[0][1][1]));
	int i = 0;
	while (!(lev == -1 && tag[0] == 1)){//See pp.148 for algo

		while(branchTerm(poldP, lev, tag, endpt, fixRep, word, imgArr) == 0){
			goForward(plev, tag, word, gens);	
			printWord(lev, tag);
		}
		do{
			goBackwards(plev);
			printWord(lev, tag);
		}while(!(availableTurn(plev, tag) == 1 || lev == -1));

		//Might need to put another exit condition check right here...
		//printf("LEV = %d\n", lev);
		if (lev == -1 && tag[0] == 1) break;	
		turnForward(plev, tag, word, gens);
		printWord(lev, tag);
	}
	saveArrayAsBMP(imgArr,filename, WIDTH, HEIGHT);

}

void main(){
	double complex ta = 2.2;
	double complex tb = 2.2;
	int numIm = 0;
	float theta = 0;

	srand(time(NULL));

	float *** imgArr = (float***)malloc(WIDTH*sizeof(float**));
	for (int i = 0; i< WIDTH; i++) {
		imgArr[i] = (float **) malloc(HEIGHT*sizeof(float *));
		for (int j = 0; j < HEIGHT; j++) 
			imgArr[i][j] = (float *) malloc(3 *sizeof(float));
	}

	while(1){

		for (int i = 0; i< WIDTH; i++) {
			for (int j = 0; j < HEIGHT; j++) {
				imgArr[i][j][0] = 0;
				imgArr[i][j][1] = 0;
				imgArr[i][j][2] = 0;
			}
		}



		//ta = 2 + 2 * I * -sin(theta);
		//tb = 2*sin(theta) + I * cos(theta);
		//ta = 1.91 + cos(theta)*I;
		//tb = 1.91 + sin(theta)*I;
		ta  = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
		tb  = (float)rand()/(float)(RAND_MAX/2) +  -I +(float)rand()/(float)(RAND_MAX/2)*I;
		//ta = 2 + cos(theta) + 2*I + I*sin(theta);
		//tb = 2 + cos(theta) + 2*I + I*sin(theta);
		computeDepthFirst(ta, tb, imgArr, numIm);
		numIm++;
		//	LEVMAX++;
		theta += 0.05;
		if (theta > 3.1415928) exit(1);
	}
}
