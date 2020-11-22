#include <stdio.h>
#include <complex.h>

#include "include/accidents.h"


double complex traceRecursion(int p, int q, double complex ta, double complex tB, double complex taB){
	//Implementation of the trace recursion algo from pp. 286
	//Here, pN and qN designs a specific numerator and denominator in a rational number-> p/q	
	//Basically, computes a polynom specific to a given p/q rational in terms of traces ta tB taB 
	if (p == 0 && q == 1) return ta;
	else if (p == 1 && q == 0) return tB;
		
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

	while(p3*q != p*q3){
		if (p * q3 < p3 * q){
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
double complex findMu(int denom){
	//This 
}

double newtPQ(int p, int q, int denom){

}

double newtonSolver(double complex initGuess){

}
