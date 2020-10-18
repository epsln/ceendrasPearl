#ifndef arraysOps
#define arraysOps
void matrix3dto3D(double complex matrix3d[1000][2][2], double complex destination[1000][2][2], int i, int k);
void matrix3dto2D(double complex matrix3d[1000][2][2], double complex destination[2][2], int i);
void matrix2dto3D(double complex matrix2d[2][2], double complex destination[1000][2][2], int i);

#endif
