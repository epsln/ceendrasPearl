#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <time.h>

#include "include/arraysOps.h"
#include "include/plot.h"

double map(double n,double  start1,double  stop1,double  start2,double  stop2){//map a real from one range to another
	return ((n-start1)/(stop1-start1))*(stop2-start2)+start2;
}

double complex randomComplex(double complex min, double complex max){
	double realPart = (creal(max - min)) * ((((double) rand()) / (double) RAND_MAX)) + creal(min) ;
	double imagPart = (cimag(max - min)) * ((((double) rand()) / (double) RAND_MAX)) + cimag(min) ;

	return realPart + I * imagPart;
}

double complex mobiusOnPoint(double complex T[2][2], double complex z){//See pp.75
	return (T[0][0] * z + T[1][0])/(T[0][1] * z + T[1][1]);
}

int modulo(int a, int b){
	return ((a % b) + b ) % b;
}

double computeBoxdim(image_t *img){
	long int count = 0;
	double epsi = 2.0 * img->bounds/(float)img->w * img->w/(float)img->h;//Might need to take into account other viewpoints...
	if (img->bitwise == 1){
		printf("Bitwise boxdim is not yet implemeted !\nExiting...\n");
		exit(-2);
	}
	else{
		for(int i = 0; i < img->w; i++){
			for(int j = 0; j < img->h; j++){
				if (img->pointArr[j * img->h + i] == 1) count++;
			}
		}
	}

        memset(img->pointArr, 0, (img->w*img->h) * (sizeof *img->pointArr));
	return log(count)/log(1/epsi);
}

void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]){//Not implementing any higher dims lol
	//(a e + b g | a f + b h
	// c e + d g | c f + d h)
	C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[0][1];
	C[1][0] = A[0][0] * B[1][0] + A[1][0] * B[1][1];
	C[0][1] = A[0][1] * B[0][0] + A[1][1] * B[0][1];
	C[1][1] = A[0][1] * B[1][0] + A[1][1] * B[1][1];
}

double complex fix(double complex T[2][2]){//See pp.76
	double complex z0 = (T[0][0] - T[1][1] - csqrt(cpow(T[1][1] - T[0][0], 2) + 4*T[1][0]*T[0][1]))/(2*T[0][1]);
	return z0;
}


void computeRepetends(double complex* gens, double complex fixRep[4][3]){//See pp.218
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



void computeRepetendsv2(double complex* gens, double complex fixRep[4][4]){//See pp.218
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

	//aBa
	matmul(buff_gen_a, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][1] = fix(buff_out1);

	//Baa
	matmul(buff_gen_B, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][2] = fix(buff_out1);

	//BAba
	matmul(buff_gen_B, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][3] = fix(buff_out0);

	//ABab
	matmul(buff_gen_A, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][0] = fix(buff_out0);

	//AAb
	matmul(buff_gen_A, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	fixRep[1][1] = fix(buff_out1);

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

	//AbA
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][1] = fix(buff_out1);

	//bAA
	matmul(buff_gen_b, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][2] = fix(buff_out1);

	//baBA
	matmul(buff_gen_b, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][3] = fix(buff_out0);

	//abAB
	matmul(buff_gen_a, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][0] = fix(buff_out0);

	//aaB
	matmul(buff_gen_a, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	fixRep[3][1] = fix(buff_out1);

	//AbaB
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][2] = fix(buff_out0);

}

void computeCycles(double complex begpt[4], double complex endpt[4], double complex* gens){

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

	//Compute the fix points of all right most turns 
	matmul(buff_gen1, buff_gen2, buff_out0); 
	matmul(buff_out0, buff_gen3, buff_out1); 
	matmul(buff_out1, buff_gen0, buff_out0); 
	begpt[0] = fix(buff_out0);

	matmul(buff_gen2, buff_gen3, buff_out0); 
	matmul(buff_out0, buff_gen0, buff_out1); 
	matmul(buff_out1, buff_gen1, buff_out0); 
	begpt[1] = fix(buff_out0);

	matmul(buff_gen3, buff_gen0, buff_out0); 
	matmul(buff_out0, buff_gen1, buff_out1); 
	matmul(buff_out1, buff_gen2, buff_out0); 
	begpt[2] = fix(buff_out0);

	matmul(buff_gen0, buff_gen1, buff_out0); 
	matmul(buff_out0, buff_gen2, buff_out1); 
	matmul(buff_out1, buff_gen3, buff_out0); 
	begpt[3] = fix(buff_out0);

}

