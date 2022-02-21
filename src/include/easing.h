#ifndef easing 
#define easing 

double easeInOutQuad(double t, double b, double c, double d);
double schlickEase(double x, double s, double t, double beg, double end, double nsteps);
double complex InOutQuadComplex(double t, double complex beg, double complex end, double nsteps);
double complex schlickComplex(double x, double s, double t, double complex beg, double complex end, double nsteps);

double complex bezier(double complex p0, double complex p1, double complex p2, double complex p3, float t);

#endif 
