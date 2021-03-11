#ifndef quatsMath
#define quatsMath 

typedef struct quat_t{
	long double a, b, c, d;
}quat_t;

void quatAdd(quat_t* out, quat_t a, quat_t b);
void quatSub(quat_t* out, quat_t a, quat_t b);
void quatMult(quat_t* out, quat_t a, quat_t b);
void quatDiv(quat_t* out, quat_t a, quat_t b);

long double quatNorm(quat_t h);
#endif
