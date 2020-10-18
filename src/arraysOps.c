#include <complex.h>
void matrix3dto3D(double complex matrix3d[1000][2][2], double complex destination[1000][2][2], int i, int k){
	destination[i][0][0] = matrix3d[k][0][0];
	destination[i][0][1] = matrix3d[k][0][1];
	destination[i][1][0] = matrix3d[k][1][0];
	destination[i][1][1] = matrix3d[k][1][1];
}

void matrix3dto2D(double complex matrix3d[1000][2][2], double complex destination[2][2], int i){
	destination[0][0] = matrix3d[i][0][0];
	destination[0][1] = matrix3d[i][0][1];
	destination[1][0] = matrix3d[i][1][0];
	destination[1][1] = matrix3d[i][1][1];
}

void matrix2dto3D(double complex matrix2d[2][2], double complex destination[1000][2][2], int i){
	destination[i][0][0] = matrix2d[0][0];
	destination[i][0][1] = matrix2d[0][1];
	destination[i][1][0] = matrix2d[1][0];
	destination[i][1][1] = matrix2d[1][1];
}


