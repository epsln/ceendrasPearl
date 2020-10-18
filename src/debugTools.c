#include <stdio.h>
#include <complex.h>
void showMatrix(double complex mat[2][2], double* PARAMS){
	int DEBUG = (int)PARAMS[6];
	if (DEBUG == 1)
		printf("[[%lf + i %lf, %lf + i%lf],\n[%lf + i %lf, %lf + i %lf ]]\n", creal(mat[0][0]),cimag(mat[0][0]), creal(mat[1][0]), cimag(mat[1][0]), creal(mat[0][1]),cimag(mat[0][1]), creal(mat[1][1]),cimag(mat[1][1]));

}


void printWord(int lev, int* tag, double* PARAMS){
	int DEBUG = (int)PARAMS[6];
	if (DEBUG == 1){
		for (int i = 0; i <= lev; i++){
			if (tag[i] == 0)
				printf("a");
			else if ( tag[i] == 1)
				printf("b");
			else if (tag[i] == 2)
				printf("A");
			else if (tag[i] == 3)
				printf("B");
			else
				printf("%d", tag[i]);
		}
		printf("\n");
	}
}


