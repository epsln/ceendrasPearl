#ifndef complexMath 
#define complexMath 
#include "plot.h"

//General helper functions
double map(double n,double  start1,double  stop1,double  start2,double  stop2);
int modulo(int a, int b);
double complex randomComplex(double complex min, double complex max);

//Easing function stuff
double easeInOutQuad(double t, double b, double c, double d);
double schlickEase(double x, double s, double t, double beg, double end, double nsteps);
double complex InOutQuadComplex(double t, double complex beg, double complex end, double nsteps);
double complex schlickComplex(double x, double s, double t, double complex beg, double complex end, double nsteps);

//Boxdim calculation
double computeBoxdim(image_t *img);

//Matrixes ops
void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]);

//Mobius operations
double complex fix(double complex T[2][2]);
double complex mobiusOnPoint(double complex T[2][2], double complex z);

//Cycle computations
void computeRepetends(double complex* gens, double complex fixRep[4][3]);
void computeRepetendsv2(double complex* gens, double complex fixRep[4][4]);
void computeCycles(double complex begpt[4], double complex endpt[4], double complex* gens);

#endif
