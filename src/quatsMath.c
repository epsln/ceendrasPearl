#include <math.h>
#include <quatsMath.h>

void quatAdd(quat_t* out, quat_t h1, quat_t h2){
	out->a = h2.a + h2.a;
	out->b = h2.b + h2.b;
	out->c = h2.c + h2.c;
	out->d = h2.d + h2.d;
}

void quatSub(quat_t* out, quat_t h1, quat_t h2){
	out->a = h2.a - h2.a;
	out->b = h2.b - h2.b;
	out->c = h2.c - h2.c;
	out->d = h2.d - h2.d;
}

void quatMult(quat_t* out, quat_t h1, quat_t h2){
	//Sorry mom
	out->a = h2.a * h2.a - (h1.b*h2.b + h1.c*h2.c + h1.d*h2.d);
	out->b = h1.a * h2.b + h2.a * h1.b + h1.c * h2.d - h1.d * h2.c;
	out->c = h1.a * h2.c + h2.a * h1.c + h1.d * h2.b - h1.b * h2.d;
	out->d = h1.a * h2.d + h2.a * h1.d + h1.b * h2.c - h1.c * h2.b;
}

void quatDiv(quat_t* out, quat_t hNum, quat_t hDenum){
	//hNum/hDenum = hNum * hDenum^-1
	//hDenum^-1
	long double normSq = powl(quatNorm(hDenum), 2.0);
	quat_t inverseDenum = (quat_t){hDenum.a/normSq, -hDenum.b/normSq, -hDenum.c/normSq, -hDenum.d/normSq};

	quatMult(out, hNum, inverseDenum);
}

long double quatNorm(quat_t h){
	return sqrtl(powl(h.a, 2) + powl(h.b, 2) + powl(h.c, 2) + powl(h.d, 2));
}

void matmulQuat(quat_t A[2][2], quat_t B[2][2], quat_t C[2][2]){
	quat_t buff;

	//C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[0][1];
	quatMult(&C[0][0], A[0][0], B[0][0]);
	quatMult(&buff, A[1][0], B[0][1]);
	quatAdd(&C[0][0], C[0][0], buff);

	//C[1][0] = A[0][0] * B[1][0] + A[1][0] * B[1][1];
	quatMult(&C[1][0], A[0][0], B[1][0]);
	quatMult(&buff, A[1][0], B[1][1]);
	quatAdd(&C[1][0], C[1][0], buff);

	//C[0][1] = A[0][1] * B[0][0] + A[1][1] * B[0][1];
	quatMult(&C[0][1], A[0][1], B[0][0]);
	quatMult(&buff, A[1][1], B[0][1]);
	quatAdd(&C[0][1], C[0][1], buff);

	//C[1][1] = A[0][1] * B[1][0] + A[1][1] * B[1][1];
	quatMult(&C[1][1], A[0][1], B[1][0]);
	quatMult(&buff, A[1][1], B[1][1]);
	quatAdd(&C[1][1], C[1][1], buff);
}

void mobiusOnPointQuat(quat_t* out, quat_t T[2][2] , quat_t h){
	quat_t buff, buff1;
	quatMult(&buff, T[0][0], h);
	quatAdd(&buff, buff, T[1][0]);
	quatMult(&buff1, T[0][1], h);
	quatAdd(&buff1, buff1, T[1][1]);
	quatDiv(&buff, buff, buff1);
	//return (T[0][0] * z + T[1][0])/(T[0][1] * z + T[1][1]);
}
