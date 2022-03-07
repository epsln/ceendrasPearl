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
	
	/*
	double controlPoints[2][4];
	controlPoints[0][0] = creal(p0); controlPoints[1][0] = cimag(p0);  
	controlPoints[0][1] = creal(p1); controlPoints[1][1] = cimag(p1);  
	controlPoints[0][2] = creal(p2); controlPoints[1][2] = cimag(p2);  
	controlPoints[0][3] = creal(p3); controlPoints[1][3] = cimag(p3);  

	px = pow(1 - t, 3) * controlPoints[0][0];
	py = pow(1 - t, 3) * controlPoints[1][0];

	px += 3 * pow(1 - t, 2) * t * controlPoints[0][1];
	py += 3 * pow(1 - t, 2) * t * controlPoints[1][1];
	
	px += 3 * (1 - t) * pow(t, 2) * controlPoints[0][1];
	py += 3 * (1 - t) * pow(t, 2) * controlPoints[1][2];

	px += pow(t, 3) * controlPoints[0][3];
	py += pow(t, 3) * controlPoints[1][3];
	*/
	printf("%d\n", n_steps);
	for (int i = 0; i < n_steps; i++){
		px += binomial(n_steps, i) * pow(1 - t, n_steps - i) * pow(t, i) * creal(controlPoints[i]); 
		py += binomial(n_steps, i) * pow(1 - t, n_steps - i) * pow(t, i) * cimag(controlPoints[i]); 
		printf("t: %lf, %lf + %lf\n", t, px, py);
		printf("(%d %d) = %d\n", n_steps, i, binomial(n_steps, i));
		printf("P_%d  = %lf %lf\n", i, creal(controlPoints[i]), cimag(controlPoints[i]));
	}
	printf("%f + %f\n", px, py);
	return px + I * py;
}
