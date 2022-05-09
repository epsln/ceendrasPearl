#include <complex.h>
#include <math.h>
#include <stdio.h>
#include "include/complexMath.h"

//Easing functions 
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

double complex bezier(float t, int n_steps, double complex controlPoints[n_steps]){
	double px = 0;
	double py = 0;
	
	int n = n_steps - 1;
	for (int i = 0; i < n_steps; i++){
		px += binomial(n, i) * pow(1 - t, n - i) * pow(t, i) * creal(controlPoints[i]); 
		py += binomial(n, i) * pow(1 - t, n - i) * pow(t, i) * cimag(controlPoints[i]); 
	}
	return px + I * py;
}
