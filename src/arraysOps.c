#include <complex.h>
#include "mpc.h"
void matrix3dto3D(mpc_t* matrix3d, mpc_t* destination, int i, int k){
	mpc_set(destination[(i * 2 + 0) * 2 + 0], matrix3d[(k * 2 + 0) * 2 + 0], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 1) * 2 + 0], matrix3d[(k * 2 + 1) * 2 + 0], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 0) * 2 + 1], matrix3d[(k * 2 + 0) * 2 + 1], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 1) * 2 + 1], matrix3d[(k * 2 + 1) * 2 + 1], MPC_RNDNN);
}

void matrix3dto2D(mpc_t* matrix3d, mpc_t destination[2][2], int i){
	mpc_set(destination[0][0], matrix3d[(i * 2 + 0) * 2 + 0], MPC_RNDNN);
	mpc_set(destination[0][1], matrix3d[(i * 2 + 0) * 2 + 1], MPC_RNDNN);
	mpc_set(destination[1][0], matrix3d[(i * 2 + 1) * 2 + 0], MPC_RNDNN);
	mpc_set(destination[1][1], matrix3d[(i * 2 + 1) * 2 + 1], MPC_RNDNN);
}

void matrix2dto3D(mpc_t matrix2d[2][2], mpc_t* destination, int i){
	mpc_set(destination[(i * 2 + 0) * 2 + 0], matrix2d[0][0], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 0) * 2 + 1], matrix2d[0][1], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 1) * 2 + 0], matrix2d[1][0], MPC_RNDNN);
	mpc_set(destination[(i * 2 + 1) * 2 + 1], matrix2d[1][1], MPC_RNDNN);
}


