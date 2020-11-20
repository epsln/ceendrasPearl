#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>

#include "include/arraysOps.h"

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

void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]){//Not implementing any higher dims lol

		//(a e + b g | a f + b h
	// c e + d g | c f + d h)
	// a = A[0][0] b = A[1][0] c = A[0][1] d = A[1][1]
	// e = B[0][0] f = B[1][0] g = B[0][1] h = B[1][1]

	C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[0][1];
	C[1][0] = A[0][0] * B[1][0] + A[1][0] * B[1][1];
	C[0][1] = A[0][1] * B[0][0] + A[1][1] * B[0][1];
	C[1][1] = A[0][1] * B[1][0] + A[1][1] * B[1][1];
}

double complex fix(double complex T[2][2]){//See pp.76
	// a = A[0][0] b = A[1][0] c = A[0][1] d = A[1][1]
	double complex z0 = (T[0][0] - T[1][1] - csqrt(cpow(T[1][1] - T[0][0], 2) + 4*T[1][0]*T[0][1]))/(2*T[0][1]);
	return z0;
}


//double easeInOutQuad (double t, double b, double c, int d) { // Page 211
//    if ((t/=d/2) < 1) return c/2*t*t + b;
//    return -c/2 * ((t - 1)*(t-2) - 1) + b;
//};


double easeInOutQuad(double t, double b, double c, double d) {
	t /= d/2;
	if (t < 1) return c/2*t*t + b;
	t--;
	return -c/2 * (t*(t-2) - 1) + b;
}

double schlickEase(double x, double s, double t, double beg, double end, double nsteps){
	//Trying to implement this paper for shits and giggles: https://arxiv.org/pdf/2010.09714.pdf
	x /= nsteps;
	if (x < t){
		double res = t*x/(x + s*(t - x) + 0.00001); 
		//printf("x: %f t:%f s:%f res: %f\n", x, t, s, res);
		return map(res, 0, 1, beg, end); 
	}
	else{
		double res = ((1 - t) * (x - 1))/(1 - x - s*(t - x) + 0.00001); 
		//printf("x: %f t:%f s:%f res: %f\n", x, t, s, res);
		return map(res, 0, 1, beg, end); 
	}
}

double complex InOutQuadComplex(double t, double complex beg, double complex end, double nsteps){
	return easeInOutQuad(t, creal(beg), creal(end), nsteps) + I * easeInOutQuad(t, cimag(beg), cimag(end), nsteps);
}

double complex schlickComplex(double x, double s, double t, double complex beg, double complex end, double nsteps){
	return schlickEase(x, s, t, creal(beg), creal(end), nsteps) + I * schlickEase(x, s, t, cimag(beg), cimag(end), nsteps);
}
void maskitRecipe(double complex ta, double complex gens[4][2][2]){
	//See pp. 259
		
	gens[0][0][0] = ta;
	gens[0][1][0] = -I; 
	gens[0][0][1] = -I;
	gens[0][1][1] = 0; 

	gens[1][0][0] = 1;
	gens[1][1][0] = 2; 
	gens[1][0][1] = 0;
	gens[1][1][1] = 1; 

	gens[2][0][0] = gens[0][1][1];
	gens[2][1][0] = -gens[0][1][0];
	gens[2][0][1] = -gens[0][0][1];
	gens[2][1][1] = gens[0][0][0];

	gens[3][0][0] = gens[1][1][1];
	gens[3][1][0] = -gens[1][1][0];
	gens[3][0][1] = -gens[1][0][1];
	gens[3][1][1] = gens[1][0][0];
}


void grandmaRecipe(double complex ta, double complex tb, double complex gens[4][2][2]){
	double complex a = 1;
	double complex b = (-ta * tb);
	double complex c = ta * ta + tb * tb; double complex delta = b*b - 4 * a * c; 
	double complex tab = (- b - csqrt(delta))/(2 * a); 
	double complex z0 = ((tab - 2) * tb)/(tb * tab - 2 * ta + 2 * I * tab);

	gens[0][0][0] = ta/2;
	gens[0][1][0] =  (ta*tab - 2 * tb + 4 * I)/(z0*(2 * tab + 4)); 
	gens[0][0][1] = ((ta * tab - 2 * tb - 4 * I)*z0)/(2* tab - 4);
	gens[0][1][1] = ta/2; 

	gens[1][0][0] = (tb - 2 * I)/2;
	gens[1][1][0] = tb/2; 
	gens[1][0][1] = tb/2;
	gens[1][1][1] = (tb + 2 * I)/2; 

	gens[2][0][0] = gens[0][1][1];
	gens[2][1][0] = -gens[0][1][0];
	gens[2][0][1] = -gens[0][0][1];
	gens[2][1][1] = gens[0][0][0];

	gens[3][0][0] = gens[1][1][1];
	gens[3][1][0] = -gens[1][1][0];
	gens[3][0][1] = -gens[1][0][1];
	gens[3][1][1] = gens[1][0][0];
}

void grandmaSpecialRecipe(double complex ta, double complex tb, double complex tab, double complex gens[4][2][2]){
	//Grandma's Special four-alarm two generator group recipe 
	//See pp. 261
	double complex tc = ta * ta + tb * tb + tab * tab - ta * tb * tab - 2;
	double complex Q  = csqrt(2 - tc);	
	double complex R = 0;
	if (cabs(tc + I * Q * csqrt(tc + 2)) >= 2)
		R = csqrt(tc + 2);
	else
		R = -csqrt(tc + 2);
	
	double complex z0 = ((tab - 2) * (tb + R))/(tb * tab - 2 * ta + I * Q * tab);
	
	gens[0][0][0] = ta/2;
	gens[0][1][0] =  (ta * tab - 2 * tb + 2 * I * Q)/(z0 * (2 * tab + 4)); 
	gens[0][0][1] = ((ta * tab - 2 * tb - 2 * I * Q) * z0)/(2 * tab - 4);
	gens[0][1][1] = ta/2; 

	gens[1][0][0] = (tb - I * Q)/2;
	gens[1][1][0] = (tb * tab - 2 * ta + I * Q * tab)/(z0*(2 * tab + 4));
	gens[1][0][1] = ((ta * tab - 2 * ta - I * Q * tab)*z0)/(2 * tab - 4);
	gens[1][1][1] = (tb + 2 * I * Q)/2; 


	gens[2][0][0] = gens[0][1][1];
	gens[2][1][0] = -gens[0][1][0];
	gens[2][0][1] = -gens[0][0][1];
	gens[2][1][1] = gens[0][0][0];

	gens[3][0][0] = gens[1][1][1];
	gens[3][1][0] = -gens[1][1][0];
	gens[3][0][1] = -gens[1][0][1];
	gens[3][1][1] = gens[1][0][0];

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



void computeRepetendsv2(double complex gens[4][2][2], double complex fixRep[4][4]){//See pp.218
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

void computeCycles(double complex begpt[4], double complex endpt[4], double complex gens[4][2][2]){

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

