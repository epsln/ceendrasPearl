#ifndef quatsMath
#define quatsMath 

typedef struct quat_t{
	long double a, b, c, d;
}quat_t;


quat_t quatAdd(quat_t a, quat_t b);
quat_t quatSub(quat_t a, quat_t b);
quat_t quatMult(quat_t a, quat_t b);
quat_t quatDiv(quat_t a, quat_t b);

long double quatNorm(quat_t h);
long double quatDist(quat_t h1, quat_t h2);

void matmulQuat(quat_t A[2][2], quat_t B[2][2], quat_t C[2][2]);
quat_t mobiusOnPointQuat(quat_t T[2][2] , quat_t h);

quat_t fixQuat(quat_t T[2][2]);
void computeRepetendsQuats(quat_t* gens, quat_t fixRep[4][4]);

void showMatrixQuat(quat_t T[2][2]);
#endif
