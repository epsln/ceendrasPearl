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
