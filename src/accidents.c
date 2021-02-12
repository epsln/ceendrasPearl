#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "include/accidents.h"


double complex tracePoly(ratio fraction, double complex ta, double complex tB, double complex taB){
	//Implementation of the trace recursion algo from pp. 286
	//Here, pN and qN designs a specific numerator and denominator in a rational number-> p/q	
	//Basically, computes a polynom specific to a given p/q rational in terms of traces ta tB taB 
	if (fraction.p == 0 && fraction.q == 1) return ta;
	else if (fraction.p== 1 && fraction.q== 0){
		return tB;
	}

		
	int p1, p2, p3;
	int q1, q2, q3;

	double complex tr_u, tr_v, tr_uv, temp;

	p1 = 0;
	q1 = 1;

	p2 = 1;
	q2 = 0;

	p3 = 1;
	q3 = 1;

	tr_u  = ta;
	tr_v  = tB;
	tr_uv = taB;

	while(p3*fraction.q!= fraction.p*q3){
		if (fraction.p* q3 < p3 * fraction.q){
			p2 = p3;
			q2 = q3;
			p3 = (p1 + p3);
			q3 = (q1 + q3);

			temp = tr_uv;

			tr_uv = tr_u * tr_uv - tr_v;
			tr_v = temp;
		}
		else{
			p1 = p3;
			q1 = q3;
			p3 = (p2 + p3);
			q3 = (q2 + q3);

			temp = tr_uv;

			tr_uv = tr_v * tr_uv - tr_u;
			tr_u = temp;
		}
	}
	return tr_uv;
}


void makeFareySeq(int denum, ratio* fareyArr){
	//Code for RosettaCode: https://rosettacode.org/wiki/Farey_sequence#C, with minor modifications to put it all in an array
	ratio f1 = {0, 1};
	ratio f2 = {1, denum};
	ratio t;
	int k = 0;
	int i = 1;

	for (int i = 0; i < denum*denum; i++){
		fareyArr[i] = (ratio) {0, 0};	
	}

	fareyArr[0] = (ratio) {0, 1};
	fareyArr[1] = (ratio) {1, denum};
	while (f2.q > 1){
		i++;
		k = (denum + f1.q)/f2.q;
		t = f1;
		f1 = f2;
		f2 = (ratio) {f2.p * k - t.p, f2.q * k - t.q};
		fareyArr[i] = f2;
	}
}

void makeFiboSeq(int lengthAnim, ratio* fiboFracts){
	long long int fiboSerie[70] = {0};//Have to stop at 70 because of long long int MAXINT being kinda small in comparison to \infty
	fiboSerie[0] = 0;
	fiboSerie[1] = 1;
	for (int i = 2; i < 69; i++){
		fiboSerie[i] = fiboSerie[i - 2] + fiboSerie[i - 1]; 
		fiboFracts[i - 2] = (ratio){fiboSerie[i - 2], fiboSerie[i]};  	
	}
}

void makePiSeq(int lengthAnim, ratio* piFracts){
	piFracts[0] = (ratio){4, 1};
	for (int i = 1; i < lengthAnim ; i++){
		piFracts[i] = (ratio){piFracts[i - 1].p * (2 * i + 1) + pow(-1, i) * piFracts[i - 1].q * 4,
	       		              piFracts[i - 1].q * (2 * i + 1)};
	}
	
}

void makeContinuedFraction(int lengthArr, double realNum, ratio* fractionArr){
	//Creates an array filed with the continued fraction approximation of realNum
	//TODO: Make this using arbitrary precision...
	long double intPart;
	long double reciprocal, num, denum, temp;
	int fractTerms[lengthArr + 1];//Array to hold the values of the denominator of the continued fraction


	intPart = floor(realNum);
	reciprocal = 1.0/(realNum - intPart);
	printf("rec: %llf\n", reciprocal);
	fractTerms[0] = intPart;		
	for (int i = 1; i < lengthArr; i++){
		intPart = floor(reciprocal);
		reciprocal = 1.0/(reciprocal - intPart);
		fractTerms[i] = (int)intPart;
		printf("rec: %llf\n", reciprocal);
		//printf("%d ", fractTerms[i]);
	}
	printf("\n");

	num = 1;
	denum = fractTerms[lengthArr - 1];
	for(int i = lengthArr - 1; i > 0; i--){
		temp = denum;
		denum = fractTerms[i] * denum + num;
		num = temp;
		fractionArr[lengthArr - i - 1] = (ratio){num, denum};	
		printf("%lld/%lld\n", fractionArr[lengthArr - i - 1].p, fractionArr[lengthArr - i - 1].q);
	}
	

}

void nextPQ(int* pP, int* pQ, int denom){
	//Doesnt work ?? ;_;
	float p1, p2, r;
	float q1, q2, s;
	int sign = -1;
	float a, k, temp0, temp1;
	p1 = 0;
	q1 = 1;
	p2 = 1;
	q2 = 0;
	r = *pP;
	s = *pQ;
	//printf("p/q: %d/%d\ndenom: %d\n", *pP, *pQ, denom);
	while (s != 0){
		a = floor((float)r/s);
		temp0 = s;
		s = r - a*s;  
		r = temp0;
		temp0 = p2;
		temp1 = q2;
		p2 = a*p2 + p1;
		q2 = a*q2 + q1;
		p1 = temp0;
		q1 = temp1;
		sign = -sign;
	//	printf("s: %.1f r: %.1f p1/q1: %.1f/%.1f p2/q2: %.1f/%.1f sign: %d\n", s, r, p1, q1, p2, q2, sign);
	}
	k = floor(((float)denom - sign*q1)/denom);
//	printf("k: %d\n", k);
	*pP = k * (*pP) + sign * p1;
	*pQ = k * (*pQ) + sign * q1; 

}

double complex traceEqn(ratio fraction, double complex mu){
	//Helper function to reduce a bit the line size :)
	//printf("tPoly(%d, %d, %lf + i %lf, %lf, %lf + I %lf) = ", fraction.p, fraction.q, creal(-I*mu), cimag(-I*mu), 2.0, creal(-I * mu + 2 * I), cimag(-I * mu + 2 * I));
	//printf("%lf + %lf\n", creal(tracePoly(fraction, -I*mu, 2, -I * mu + 2 * I) - 2), cimag(tracePoly(fraction, -I*mu, 2, -I * mu + 2 * I) - 2));
	return (tracePoly(fraction, -I*mu, 2, -I * mu + 2 * I) - 2);
}

void newtonSolver(double complex *pz0, ratio fraction){
	//TODO: Maybe implement the halley ? 
	int maxiter = 10000;
	double complex z = *pz0;
	double epsilon = 1E-15;
	double complex realVal, imagVal, deriv;
	double complex traceEqVal;
	
	//Carefull, without a imaginary part, newton doesnt converge !
	*pz0 += I;
	for (int i = 0; i < maxiter; i++){
		//Compute the complex derivate using a simple finite differences scheme on real and imag axis
		realVal = (traceEqn(fraction, z + epsilon)     - traceEqn(fraction, z - epsilon))/(2*epsilon);
		imagVal = (traceEqn(fraction, z + epsilon * I) - traceEqn(fraction, z - epsilon * I))/(2 * epsilon * I);
		deriv = (realVal + imagVal)/2.0;
		//Update the guess 	
		traceEqVal = traceEqn(fraction, z);
		z = z - traceEqVal/deriv;
		if (cabs(traceEqVal) <= 1E-7 && cabs(z - *pz0) <= 1E-7){
			return;
		}
		else
			*pz0 = z;
	}
	printf("Newton method failed !\nAbandoning all hopes and exiting...\n");
	exit(-1);
}
