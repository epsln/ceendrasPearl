#ifndef arraysOps
#define arraysOps
void matrix3dto3D(mpc_t* matrix3d, mpc_t *destination, int i, int k);
void matrix3dto2D(mpc_t* matrix3d, mpc_t destination[2][2], int i);
void matrix2dto3D(mpc_t matrix2d[2][2], mpc_t* destination, int i);

#endif
