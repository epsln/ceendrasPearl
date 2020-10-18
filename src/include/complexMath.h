#ifndef complexMath 
#define complexMath 

double complex mobiusOnPoint(double complex T[2][2], double complex z);
int modulo(int a, int b);
void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]);
double complex fix(double complex T[2][2]);
void grandmaRecipe(double complex ta, double complex tb, double complex gens[4][2][2]);
void computeRepetends(double complex gens[4][2][2], double complex fixRep[4][3]);

#endif
