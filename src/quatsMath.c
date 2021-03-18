#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <quatsMath.h>
#include <arraysOps.h>


quat_t quatAdd(quat_t h1, quat_t h2){
	return (quat_t){h1.a + h2.a,
			h1.b + h2.b,
			h1.c + h2.c,
			h1.d + h2.d
	};
}

quat_t quatSub(quat_t h1, quat_t h2){
	return (quat_t){h1.a - h2.a,
			h1.b - h2.b,
			h1.c - h2.c,
			h1.d - h2.d
	};
}

quat_t quatMult(quat_t h1, quat_t h2){
	//Sorry mom
	return (quat_t){
		h2.a * h2.a - (h1.b*h2.b + h1.c*h2.c + h1.d*h2.d),
		h1.a * h2.b + h2.a * h1.b + h1.c * h2.d - h1.d * h2.c,
		h1.a * h2.c + h2.a * h1.c + h1.d * h2.b - h1.b * h2.d,
		h1.a * h2.d + h2.a * h1.d + h1.b * h2.c - h1.c * h2.b
	};
}
quat_t quatDiv(quat_t hNum, quat_t hDenum){
	//hNum/hDenum = hNum * hDenum^-1
	//hDenum^-1
	long double normSq = powl(quatNorm(hDenum), 2.0);
	quat_t inverseDenum = (quat_t){hDenum.a/normSq, -hDenum.b/normSq, -hDenum.c/normSq, -hDenum.d/normSq};

	return quatMult(hNum, inverseDenum);
}

long double quatNorm(quat_t h){
	return sqrtl(powl(h.a, 2) + powl(h.b, 2) + powl(h.c, 2) + powl(h.d, 2));
}

void matmulQuat(quat_t A[2][2], quat_t B[2][2], quat_t C[2][2]){
	quat_t buff;

	//C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[0][1];
	C[0][0] = quatMult(A[0][0], B[0][0]);
	buff    = quatMult(A[1][0], B[0][1]);
	C[0][0] = quatAdd(C[0][0], buff);

	//C[1][0] = A[0][0] * B[1][0] + A[1][0] * B[1][1];
	C[1][0] = quatMult(A[0][0], B[1][0]);
	buff    = quatMult(A[1][0], B[1][1]);
	C[1][0] = quatAdd(C[1][0], buff);

	//C[0][1] = A[0][1] * B[0][0] + A[1][1] * B[0][1];
	C[0][1] = quatMult(A[0][1], B[0][0]);
	buff    = quatMult(A[1][1], B[0][1]);
	C[0][1] = quatAdd(C[0][1], buff);

	//C[1][1] = A[0][1] * B[1][0] + A[1][1] * B[1][1];
	C[1][1] = quatMult(A[0][1], B[1][0]);
	buff    = quatMult(A[1][1], B[1][1]);
	C[1][1] = quatAdd(C[1][1], buff);
}

long double quatDist(quat_t h1, quat_t h2){
	return quatNorm(quatSub(h1, h2));
}

quat_t mobiusOnPointQuat(quat_t T[2][2], quat_t h){
	//(T[0][0] * h + T[1][0])/(T[0][1] * h + T[1][1]);
	quat_t num, denum;
	num = quatMult(T[0][0], h);
	num = quatAdd(num, T[1][0]);
	denum = quatMult(T[0][1], h);
	denum = quatAdd(denum, T[1][1]);
	return quatDiv(num, denum);
}

quat_t quatSqrt(quat_t q){
	//Thanks John D. Cook, always here when we need you
	//https://www.johndcook.com/blog/2021/01/06/quaternion-square-roots/
	long double r, theta;
	quat_t u, x;
	r = quatNorm(q);
	theta = acos(q.a/r);
	u = q;
	u.a = 0;

	x.a = pow(r, 0.5) * cos(theta/2.0);
	x.b = pow(r, 0.5) * sin(theta/2.0) * u.b;
	x.c = pow(r, 0.5) * sin(theta/2.0) * u.c;
	x.d = pow(r, 0.5) * sin(theta/2.0) * u.d;
	
	return x;	
}

quat_t fixQuat(quat_t T[2][2]){
	quat_t four = (quat_t){4, 0, 0, 0};
	quat_t two = (quat_t){2, 0, 0, 0};
	quat_t buff, buff1, num, denum;
	num = quatSub(T[0][0], T[1][1]);	
	buff = quatSub(T[1][1], T[0][0]);
	buff = quatMult(buff, buff);
	buff1 = quatMult(four, T[1][0]);
	buff1 = quatMult(buff1, T[0][1]);
	buff = quatAdd(buff, buff1);
	num = quatSub(num, buff);

	denum = quatMult(two, T[0][1]);
	
       	return quatDiv(num, denum);	
	//(T[0][0] - T[1][1] - csqrt(cpow(T[1][1] - T[0][0], 2) + 4*T[1][0]*T[0][1]))/(2*T[0][1]);
	
}

void computeRepetendsQuats(quat_t* gens, quat_t fixRep[4][4]){
	quat_t buff_gen_a[2][2];
	quat_t buff_gen_b[2][2];
	quat_t buff_gen_A[2][2];
	quat_t buff_gen_B[2][2];
	quat_t buff_out0[2][2];
	quat_t buff_out1[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2DQuat(gens, buff_gen_a, 0);
	matrix3dto2DQuat(gens, buff_gen_b, 1);
	matrix3dto2DQuat(gens, buff_gen_A, 2);
	matrix3dto2DQuat(gens, buff_gen_B, 3);

	//bAba
	matmulQuat(buff_gen_b, buff_gen_A, buff_out0);
	matmulQuat(buff_out0, buff_gen_b, buff_out1);
	matmulQuat(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][0] = fixQuat(buff_out0);

	//aBa
	matmulQuat(buff_gen_a, buff_gen_B, buff_out0);
	matmulQuat(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][1] = fixQuat(buff_out1);

	//Baa
	matmulQuat(buff_gen_B, buff_gen_a, buff_out0);
	matmulQuat(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][2] = fixQuat(buff_out1);

	//BAba
	matmulQuat(buff_gen_B, buff_gen_A, buff_out0);
	matmulQuat(buff_out0, buff_gen_b, buff_out1);
	matmulQuat(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][3] = fixQuat(buff_out0);

	//ABab
	matmulQuat(buff_gen_A, buff_gen_B, buff_out0);
	matmulQuat(buff_out0, buff_gen_a, buff_out1);
	matmulQuat(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][0] = fixQuat(buff_out0);

	//AAb
	matmulQuat(buff_gen_A, buff_gen_A, buff_out0);
	matmulQuat(buff_out0, buff_gen_b, buff_out1);
	fixRep[1][1] = fixQuat(buff_out1);

	//aBAb
	matmulQuat(buff_gen_a, buff_gen_B, buff_out0);
	matmulQuat(buff_out0, buff_gen_A, buff_out1);
	matmulQuat(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][2] = fixQuat(buff_out0);

	//BabA
	matmulQuat(buff_gen_B, buff_gen_a, buff_out0);
	matmulQuat(buff_out0, buff_gen_b, buff_out1);
	matmulQuat(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][0] = fixQuat(buff_out0);

	//AbA
	matmulQuat(buff_gen_A, buff_gen_b, buff_out0);
	matmulQuat(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][1] = fixQuat(buff_out1);

	//bAA
	matmulQuat(buff_gen_b, buff_gen_A, buff_out0);
	matmulQuat(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][2] = fixQuat(buff_out1);

	//baBA
	matmulQuat(buff_gen_b, buff_gen_a, buff_out0);
	matmulQuat(buff_out0, buff_gen_B, buff_out1);
	matmulQuat(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][3] = fixQuat(buff_out0);

	//abAB
	matmulQuat(buff_gen_a, buff_gen_b, buff_out0);
	matmulQuat(buff_out0, buff_gen_A, buff_out1);
	matmulQuat(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][0] = fixQuat(buff_out0);

	//aaB
	matmulQuat(buff_gen_a, buff_gen_a, buff_out0);
	matmulQuat(buff_out0, buff_gen_B, buff_out1);
	fixRep[3][1] = fixQuat(buff_out1);

	//AbaB
	matmulQuat(buff_gen_A, buff_gen_b, buff_out0);
	matmulQuat(buff_out0, buff_gen_a, buff_out1);
	matmulQuat(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][2] = fixQuat(buff_out0);

}

void showMatrixQuat(quat_t T[2][2]){
	printf("[%Lf %Lf %Lf %Lf],",    T[0][0].a, T[0][0].b, T[0][0].c, T[0][0].d);
	printf("[%Lf %Lf %Lf %Lf]\n",   T[1][0].a, T[1][0].b, T[1][0].c, T[1][0].d);
	printf("[%Lf %Lf %Lf %Lf],",    T[0][1].a, T[0][1].b, T[0][1].c, T[0][1].d);
	printf("[%Lf %Lf %Lf %Lf]\n\n", T[1][1].a, T[1][1].b, T[1][1].c, T[1][1].d);
}
