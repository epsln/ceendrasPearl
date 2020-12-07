#ifndef arraysOps
#define arraysOps
void matrix3dto3D(double complex* matrix3d, double complex *destination, int i, int k);
void matrix3dto2D(double complex* matrix3d, double complex destination[2][2], int i);
void matrix2dto3D(double complex matrix2d[2][2], double complex* destination, int i);

#endif
