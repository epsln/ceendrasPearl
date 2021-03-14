#ifndef arraysOps
#define arraysOps

#include "quatsMath.h"

void matrix3dto3D(double complex* matrix3d, double complex *destination, int i, int k);
void matrix3dto2D(double complex* matrix3d, double complex destination[2][2], int i);
void matrix2dto3D(double complex matrix2d[2][2], double complex* destination, int i);

void matrix3dto3DQuat(quat_t* matrix3d, quat_t *destination, int i, int k);
void matrix3dto2DQuat(quat_t* matrix3d, quat_t destination[2][2], int i);
void matrix2dto3DQuat(quat_t matrix2d[2][2], quat_t* destination, int i);

#endif
