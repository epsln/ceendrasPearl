#include <complex.h>
void matrix3dto3D(double complex* matrix3d, double complex* destination, int i, int k){
	destination[(i * 2 + 0) * 2 + 0] = matrix3d[(k * 2 + 0) * 2 + 0];
	destination[(i * 2 + 1) * 2 + 0] = matrix3d[(k * 2 + 1) * 2 + 0];
	destination[(i * 2 + 0) * 2 + 1] = matrix3d[(k * 2 + 0) * 2 + 1];
	destination[(i * 2 + 1) * 2 + 1] = matrix3d[(k * 2 + 1) * 2 + 1];
}

void matrix3dto2D(double complex* matrix3d, double complex destination[2][2], int i){
	destination[0][0] = matrix3d[(i * 2 + 0) * 2 + 0];
	destination[0][1] = matrix3d[(i * 2 + 0) * 2 + 1];
	destination[1][0] = matrix3d[(i * 2 + 1) * 2 + 0];
	destination[1][1] = matrix3d[(i * 2 + 1) * 2 + 1];
}

void matrix2dto3D(double complex matrix2d[2][2], double complex* destination, int i){
	destination[(i * 2 + 0) * 2 + 0] = matrix2d[0][0];
	destination[(i * 2 + 0) * 2 + 1] = matrix2d[0][1];
	destination[(i * 2 + 1) * 2 + 0] = matrix2d[1][0];
	destination[(i * 2 + 1) * 2 + 1] = matrix2d[1][1];
}


