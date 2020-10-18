#include <math.h>
#include <stdio.h>
#include <complex.h>

#include "include/arraysOps.h"

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
	double complex z1 = (T[0][0] - T[1][1] + csqrt(cpow(T[1][1] - T[0][0], 2) + 4*T[1][0]*T[0][1]))/(2*T[0][1]);
	return z0;
}

void grandmaRecipe(double complex ta, double complex tb, double complex gens[4][2][2]){
	printf("ta: %lf + i %lf\n", creal(ta), cimag(ta));
	printf("tb: %lf + i %lf\n", creal(tb), cimag(tb));
	double complex a = 1;
	double complex b = (-ta * tb);
	double complex c = ta * ta + tb * tb;
	double complex delta = b*b - 4 * a * c; 
	double complex tab = (- b - csqrt(delta))/(2 * a); 
	//double complex tab = ((-ta * tb) - sqrt((-ta * tb)*(-ta * tb)- 4 *  (ta * ta + tb * tb)))/(2 );
	double complex z0 = ((tab - 2) * tb)/(tb * tab - 2 * ta + 2 * I * tab);
	printf("tab = %lf + i %lf\n", creal(tab), cimag(tab));
	printf("z0 = %lf + i %lf\n", creal(z0), cimag(z0));

	double complex num = (ta*tab - 2 * tb + 4 * I);
	double complex denum = ((2 * tab + 4)*z0) ;
	printf("num = %lf + i %lf\n", creal(num), cimag(num));
	printf("denum = %lf + i %lf\n", creal(denum), cimag(denum));

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



