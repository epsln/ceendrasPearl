#ifndef readFiles 
#define readFiles 

void readPoints(float *pointsList[2]);

void readConf(int *RESX, int *RESY, int *NPOINTS, int *JITTER, float *JITTER_B, int *GRAY, int *OPENCL_ITER, int *MAXITER,int *MAXITER_G, int *MAXITER_B, float *RED_C, float *GRE_C, float *BLU_C, char filename[256], int *RAND, char kernelFilename[256]);

#endif
