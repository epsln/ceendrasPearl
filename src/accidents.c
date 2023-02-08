#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#include "include/accidents.h"
#include "include/complexMath.h"


double complex tracePoly(ratio fraction, double complex ta, double complex tB, double complex taB){
	//Implementation of the trace recursion algo from pp. 286
	//Here, pN and qN designs a specific numerator and denominator in a rational number-> p/q	
	//Basically, computes a polynom specific to a given p/q rational in terms of traces ta tB taB 
	if (fraction.p == 0 && fraction.q == 1) return ta;
	else if (fraction.p== 1 && fraction.q== 0){
		return tB;
	}


	long int p1, p2, p3;
	long int q1, q2, q3;

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

	int sizeArray =  3 * denum * denum /(3.141592653579397 * 3.141592653579397);
	for (long int idx = 0; idx < sizeArray + 10; idx++){
		fareyArr[idx] = (ratio) {0, 0};	
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
	//TODO: Rewrite me so that we don't write into an array (Just return the next element)
}

ratio getNextFareyElem(long int denum, ratio f1, ratio f2){
	ratio t;
	long int k = (denum + f1.q)/f2.q;
	t = f1;
	f1 = f2;
	return (ratio) {f2.p * k - t.p, f2.q * k - t.q};
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

ratio add_fract(ratio a, ratio b){
	return (ratio) {a.p * b.q + a.q * b.p, a.q * b.q};
}

long int gcd(long int a, long int b){
	long int temp;
	while (b != 0){
		temp = a % b;

		a = b;
		b = temp;
	}
	return a;
}

ratio simplify_fract(ratio a){
	long int divisor = gcd(a.p, a.q);	
	return (ratio) {a.p/divisor, a.q/divisor};
}

int makeContinuedFraction(int lengthArr, double realNum, ratio* fractionArr){
	//Creates an array filed with the continued fraction approximation of realNum
	double intPart;
	double reciprocal, num, denum, temp;
	long int fractTerms[lengthArr + 1];//Array to hold the values of the denominator of the continued fraction
	printf("realNum: %lf\n", realNum);
	if (realNum >= 1){
		intPart = floor(realNum);
		realNum = realNum - intPart;
	}
	intPart = floor(realNum);
	reciprocal = 1.0/(realNum - intPart);
	printf("rec: %lf int %lf %lf\n", reciprocal, intPart, realNum - intPart);
	fractTerms[0] = (long int) floor(reciprocal);
	for (int i = 1; i < lengthArr; i++){
		intPart = floor(reciprocal);
		if (reciprocal - intPart < 0.1){
			lengthArr = i;
			break;
		}
		reciprocal = 1.0/(reciprocal - intPart);
		fractTerms[i] = (long int)intPart;
		printf("fractTerms[%d] = %lld\n", i, fractTerms[i]);
	}

	if (lengthArr == 1){
		fractionArr[0] = (ratio){1, fractTerms[0]};
		printf("fractionArr[0] = %lld/%lld\n", fractionArr[0].p, fractionArr[0].q);
		return 1;
	}

	for (int i = 0; i < lengthArr; i++){
		ratio a = (ratio) {1, fractTerms[i + 1]};	
		ratio b = (ratio) {fractTerms[i], 1};	
		ratio out = add_fract(a, b);
		for (int j = i - 1; j >= 0; j--){
			if (i == 0) break;
			ratio a = (ratio) {fractTerms[j], 1};	
			ratio b = (ratio) {out.q, out.p};	
			out = add_fract(a, b);
		}
		out = (ratio) {out.q, out.p};
		fractionArr[i] = simplify_fract(out);	
		printf("fractionArr[%d] = %lld/%lld\n", i, fractionArr[i].p, fractionArr[i].q);
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

void getSpecialWordFromFract(ratio fraction, char* specialWord){
	//Fraction 2/5 -> aaaBaaB (since we don't care about the order and we'll check all cyclicPerms)
	//pp. 276
	//specialWord = calloc(fraction.p + fraction.q, sizeof(char));
	char* buff  = calloc(fraction.p + fraction.q, sizeof(char));
	int num = 1; 
	int i = 0;
	int idx = 3;
	int numToAdd = fraction.q;
	do{
		if(num + fraction.q > fraction.p + fraction.q){
			numToAdd = -fraction.p;
			idx = 0;
		}	
		if(num - fraction.p < 1){
			numToAdd = fraction.q;
			idx = 3;
		}	
		num += numToAdd;
		buff[i] = idx;
		i++;
	}while(num != 1);
	num = fraction.q + fraction.p - 1;
	for(int i = 0; i < fraction.p + fraction.q; i++){
		specialWord[i] = buff[num];
		num--;
	}
}

//void getTraceFromFract(double complex *pz0, ratio fraction){
//	//Get the trace of a particular fraction following each farey neighbor
//	int sizeArray =  3 * fraction.q * fraction.q/(3.141592653579397 * 3.141592653579397);
//	printf("sizeArray: %d\n", sizeArray);
//	ratio fractionArr[sizeArray];
//	makeFareySeq(fraction.q, fractionArr);
//	int i = -1;
//	do{
//		i++;
//		newtonSolver(pz0, fractionArr[i]);
//	}while(fraction.p != fractionArr[i].p || fraction.q != fractionArr[i].q);
//}

void getTraceFromFract(double complex *pz0, ratio fraction){
	//Get the trace of a particular fraction following each farey neighbor
	ratio f1 = (ratio) {0, 1};
	ratio f2 = (ratio) {1, fraction.q};
	ratio f2_prev = f2;
	*pz0 = 2 * I;
	newtonSolver(pz0, f1);
	while(fraction.p != f2.p || fraction.q != f2.q){
		newtonSolver(pz0, f2);
		f2 = getNextFareyElem(fraction.q, f1, f2);
		f1 = f2_prev;
		f2_prev = f2;
	}
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
	double epsilon = 1E-6;
	double complex realVal, imagVal, deriv;
	double complex traceEqVal;

	//Carefull, without a imaginary part, newton doesnt converge !
	if (cimag(*pz0) == 0)
		*pz0 += I;

	if( isinf(cimag(traceEqn(fraction, z))) || isnan(cimag(traceEqn(fraction, z)))){
		printf("Newton method failed! Getting infinity \nAbandoning all hopes and exiting...\n");
		exit(-1);
	}
	for (int i = 0; i < maxiter; i++){
		//Compute the complex derivate using a simple finite differences scheme on real and imag axis
		realVal = (traceEqn(fraction, z + epsilon)     - traceEqn(fraction, z - epsilon))/(2*epsilon);
		imagVal = (traceEqn(fraction, z + epsilon * I) - traceEqn(fraction, z - epsilon * I))/(2 * epsilon * I);
		deriv = (realVal + imagVal)/2.0;
		//Update the guess 	
		traceEqVal = traceEqn(fraction, z);
		z = z - traceEqVal/deriv;
		if (cabs(traceEqVal) <= 1E-5 && cabs(z - *pz0) <= 1E-4){
			return;
		}
		else{
			*pz0 = z;
		}
	}
	printf("z: %lf %+lf\n", creal(z), cimag(z));
	printf("Newton method failed !\nAbandoning all hopes and exiting...\n");
	exit(-1);
}
