#include <stdio.h>
#include <complex.h>

#include "include/debugTools.h"
void showMatrix(double complex mat[2][2], image_t* img){
	if (img->debug == 1)
		printf("[[%lf + i %lf, %lf + i%lf],\n[%lf + i %lf, %lf + i %lf ]]\n", creal(mat[0][0]),cimag(mat[0][0]), creal(mat[1][0]), cimag(mat[1][0]), creal(mat[0][1]),cimag(mat[0][1]), creal(mat[1][1]),cimag(mat[1][1]));

}


void printWord(int lev, int* tag, image_t* img){
	if (img->debug == 1){
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


