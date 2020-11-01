#ifndef complexMath 
#define complexMath 

//General helper functions
double map(double n,double  start1,double  stop1,double  start2,double  stop2);
int modulo(int a, int b);
double easeInOutQuad(double t, double b, double c, int d);
double complex InOutQuadComplex(double t, double complex beg, double complex end, int nsteps);

//Matrixes ops
void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]);

//Mobius operations
double complex fix(double complex T[2][2]);
double complex mobiusOnPoint(double complex T[2][2], double complex z);

//Generators computing
void grandmaRecipe(double complex ta, double complex tb, double complex gens[4][2][2]);
void grandmaSpecialRecipe(double complex ta, double complex tb, double complex tab, double complex gens[4][2][2]);

//Cycle computations
void computeRepetends(double complex gens[4][2][2], double complex fixRep[4][3]);
void computeRepetendsv2(double complex gens[4][2][2], double complex fixRep[4][4]);
void computeCycles(double complex begpt[4], double complex endpt[4], double complex gens[4][2][2]);

#endif
