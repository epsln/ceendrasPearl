#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <string.h>
#include <time.h>

#include "include/arraysOps.h"
#include "include/plot.h"

double map(double n,double  start1,double  stop1,double  start2,double  stop2){//map a real from one range to another
	return ((n-start1)/(stop1-start1))*(stop2-start2)+start2;
}

double cmap(double complex n,double complex start1, double complex stop1,double  complex start2,double complex stop2){//map a complex from one range to another
	return ((n-start1)/(stop1-start1))*(stop2-start2)+start2;
}

double complex randomComplex(double complex min, double complex max){
	double realPart = (creal(max - min)) * ((((double) rand()) / (double) RAND_MAX)) + creal(min) ;
	double imagPart = (cimag(max - min)) * ((((double) rand()) / (double) RAND_MAX)) + cimag(min) ;

	return realPart + I * imagPart;
}

double complex randomComplexFixDist(double complex z0, double d){//Get a complex number at fixed modulus and rand arg from z0
	double theta = (double)rand()/(double) RAND_MAX;
	return d * (cos(theta) + I * sin(theta)) + z0;
}

double complex mobiusOnPoint(double complex T[2][2], double complex z){//See pp.75
	return (T[0][0] * z + T[1][0])/(T[0][1] * z + T[1][1]);
}

int modulo(int a, int b){
	return ((a % b) + b ) % b;
}

double computeBoxdim(image_t *img){
	long int count = 0;
	double epsi = 2.0 * img->bounds/(float)img->w * img->w/(float)img->h;//Might need to take into account other viewpoints...
	if (img->bitwise == 1){
		printf("Bitwise boxdim is not yet implemeted !\nExiting...\n");
		exit(-2);
	}
	else{
		for(int i = 0; i < img->w; i++){
			for(int j = 0; j < img->h; j++){
				if (img->pointArr[j * img->h + i] == 1) count++;
			}
		}
	}

        memset(img->pointArr, 0, (img->w*img->h) * (sizeof *img->pointArr));
	return log(count)/log(1/epsi);
}

void matmul(double complex A[2][2], double complex B[2][2], double complex C[2][2]){//Not implementing any higher dims lol
	//(a e + b g | a f + b h
	// c e + d g | c f + d h)
	C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[0][1];
	C[1][0] = A[0][0] * B[1][0] + A[1][0] * B[1][1];
	C[0][1] = A[0][1] * B[0][0] + A[1][1] * B[0][1];
	C[1][1] = A[0][1] * B[1][0] + A[1][1] * B[1][1];
}

int checkDist(double complex* fixRep, int wordLength, double complex word[2][2], int genIndex, int numFP, double epsi){
	//Check if the distance between a word and all fixed point of the special word is < epsi
	for (int i = 0; i < numFP; i++){
		//printf("%+lf %+lf\n", creal(mobiusOnPoint(word, fixRep[genIndex * wordLength + i])), cimag(mobiusOnPoint(word, fixRep[genIndex * wordLength + i])));
		//printf("%+lf %+lf\n", creal(mobiusOnPoint(word, fixRep[genIndex * wordLength + i + 1])), cimag(mobiusOnPoint(word, fixRep[genIndex * wordLength + i + 1])));
		//printf("%lf\n",cabs(mobiusOnPoint(word, fixRep[genIndex * wordLength + i]) - mobiusOnPoint(word, fixRep[genIndex * wordLength + i + 1]))  );
		if (cabs(mobiusOnPoint(word, fixRep[genIndex * wordLength + i]) - mobiusOnPoint(word, fixRep[genIndex * wordLength + i + 1])) > epsi)
				return 0;
	}
	return 1;	
}

void composeGen(double complex* gens, int index, double complex* word){
	//Ugly workaround to those  *** arrays
	//Uses a buffer array to not have multiple buffer in calling funct
	//You get an input [2][2] mat which is also your output
	//This funct mults it by the gen at given index 
	//(a e + b g | a f + b h
	// c e + d g | c f + d h)
	// TODO: Check the indexes pls 
	double complex C[2][2]; 
	C[0][0] = word[0 * 2 + 0] * gens[(index * 2 + 0) * 2 + 0] + word[0 * 2 + 1] * gens[(index * 2 + 1) * 2 + 0];
	C[1][0] = word[0 * 2 + 0] * gens[(index * 2 + 0) * 2 + 1] + word[0 * 2 + 1] * gens[(index * 2 + 1) * 2 + 1];
	C[0][1] = word[1 * 2 + 0] * gens[(index * 2 + 0) * 2 + 0] + word[1 * 2 + 1] * gens[(index * 2 + 1) * 2 + 0];
	C[1][1] = word[1 * 2 + 0] * gens[(index * 2 + 0) * 2 + 1] + word[1 * 2 + 1] * gens[(index * 2 + 1) * 2 + 1];

	word[0 * 2 + 0] = C[0][0]; 
	word[0 * 2 + 1] = C[1][0]; 
	word[1 * 2 + 0] = C[0][1]; 	
	word[1 * 2 + 1] = C[1][1]; 
}

double complex fix(double complex T[2][2]){//See pp.76
	double complex z0 = (T[0][0] - T[1][1] - csqrt(cpow(T[1][1] - T[0][0], 2) + 4*T[1][0]*T[0][1]))/(2*T[0][1]);
	return z0;
}

void computeRepetends(double complex* gens, double complex fixRep[4][3]){//See pp.218
	double complex buff_gen_a[2][2];
	double complex buff_gen_b[2][2];
	double complex buff_gen_A[2][2];
	double complex buff_gen_B[2][2];
	double complex buff_out0[2][2];
	double complex buff_out1[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2D(gens, buff_gen_a, 0);
	matrix3dto2D(gens, buff_gen_b, 1);
	matrix3dto2D(gens, buff_gen_A, 2);
	matrix3dto2D(gens, buff_gen_B, 3);

	//bAba
	matmul(buff_gen_b, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][0] = fix(buff_out0);

	//a
	fixRep[0][1] = fix(buff_gen_a);

	//BAba
	matmul(buff_gen_B, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][2] = fix(buff_out0);

	//ABab
	matmul(buff_gen_A, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][0] = fix(buff_out0);

	//b
	fixRep[1][1] = fix(buff_gen_b);

	//aBAb
	matmul(buff_gen_a, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][2] = fix(buff_out0);

	//BabA
	matmul(buff_gen_B, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][0] = fix(buff_out0);

	//A
	fixRep[2][1] = fix(buff_gen_A);

	//baBA
	matmul(buff_gen_b, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][2] = fix(buff_out0);

	//abAB
	matmul(buff_gen_a, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][0] = fix(buff_out0);

	//B
	fixRep[3][1] = fix(buff_gen_B);

	//AbaB
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][2] = fix(buff_out0);
}

void getCyclicPerm(char** cyclicPerms, char* repr){
	//There are len(n) permutations of a string, so we store them in a 2D array size (len(n), len(n))
	for (int i = 0; i < strlen(repr); i++){
		for (int j = 0; j < strlen(repr); j++){
			cyclicPerms[i][j] = repr[j + i % strlen(repr)];
		}
	}
	
}

void makeWord(double complex out[2][2], double complex* gens, char* word, int wordLength){
	//Create a word based on its integer representation 
	double complex buff_gen_a[2][2];
	double complex buff_gen_b[2][2];
	double complex buff_gen_A[2][2];
	double complex buff_gen_B[2][2];
	double complex *buff_out = calloc(2 * 2, sizeof(double complex));

	
	//First copy out a gen to a flattened array that we'll use to compose a word
	//TODO: sort out the [2][2] and the flattened (*) array please
	//TODO: Check the indexes pls 
	

	buff_out[0 * 2 + 0] = gens[((int)word[0] * 2 + 0) * 2 + 0];
	buff_out[0 * 2 + 1] = gens[((int)word[0] * 2 + 1) * 2 + 0];
	buff_out[1 * 2 + 0] = gens[((int)word[0] * 2 + 0) * 2 + 1];
	buff_out[1 * 2 + 1] = gens[((int)word[0] * 2 + 1) * 2 + 1];

	//printf("[%lf + %lf, ",    creal(buff_out[0 * 2 + 0]), cimag(buff_out[0 * 2 + 0]));
	//printf("%lf + %lf]\n",    creal(buff_out[0 * 2 + 1]), cimag(buff_out[0 * 2 + 1]));
	//printf("[%lf + %lf,",     creal(buff_out[1 * 2 + 0]), cimag(buff_out[1 * 2 + 0]));
	//printf("[%lf + %lf]\n\n", creal(buff_out[1 * 2 + 1]), cimag(buff_out[1 * 2 + 1]));

	matrix3dto2D(gens, buff_gen_a, 0);
	matrix3dto2D(gens, buff_gen_b, 1);
	matrix3dto2D(gens, buff_gen_A, 2);
	matrix3dto2D(gens, buff_gen_B, 3);

	//Start by multiplying the 2 highest decimals
	//Length of a number in dec is log10(num) + 1
	//Decimal at index i of number is number % 
	composeGen(gens, word[0], buff_out);
	for (int i = 1; i <  wordLength; i++){
		composeGen(gens, (int)word[i] , buff_out);
	}
   
	//TODO: Check the indexes pls 
	out[0][0] = buff_out[0 * 2 + 0];
	out[1][0] = buff_out[0 * 2 + 1];
	out[0][1] = buff_out[1 * 2 + 0];
	out[1][1] = buff_out[1 * 2 + 1];
}

void computeRepetendsv2(double complex* gens, double complex* fixRep, int numFP[4], char* specialWord, int wordLength){//See pp.218
	//TODO: Clean me up
	double complex buff_gen_a[2][2];
	double complex buff_gen_b[2][2];
	double complex buff_gen_A[2][2];
	double complex buff_gen_B[2][2];
	double complex buff_out0[2][2];
	double complex buff_out1[2][2];
	double complex buffMat[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2D(gens, buff_gen_a, 0);
	matrix3dto2D(gens, buff_gen_b, 1);
	matrix3dto2D(gens, buff_gen_A, 2);
	matrix3dto2D(gens, buff_gen_B, 3);

	char* cyclicPerm = (char*) calloc(wordLength, sizeof(char));
	
	/* We have a fixRep array that is 4 * strlen(specialWord)
	 * On each index x we store all the permutations that ends in a specific generator
	 * We need to keep track of the number of permutations we add to each index of fixRep
	 * Sol: an array of index of length 4
	 * We loop on all cyclic permutations of the special words
	 * The index of the gen is the last char of the string or specialWord[len(specialWord)]
	 * We add the fixed point of this particular permutations to fixRep
	 * We increment the count of number of perms on the index corresponding to the gen once we are done
	 * We also need to account for the inverse of all special word cyclic permutations
	 * So 2 * the number of cyclic permutations
	 * If wordLength = 0 then no specialWord is defined and we just use the abAB word
	*/	
	int sizeArr = wordLength + 4;
	
	for (int i = 0; i < wordLength; i++){
		for (int j = 0; j < wordLength; j++){
			cyclicPerm[j] = specialWord[(j + i) % wordLength];
	//		printf("%d", cyclicPerm[j]);
		}
	//	printf("\n");
		makeWord(buffMat, gens, cyclicPerm, wordLength);
		fixRep[(int)cyclicPerm[wordLength - 1] * sizeArr + numFP[(int)cyclicPerm[wordLength - 1]]] = fix(buffMat);
		//printf("fix: %+lf %+lf\n", creal(fix(buffMat)), cimag(fix(buffMat)));
	//	printf("fp[%d][%d] = %lf + %lf\n",(int)cyclicPerm[wordLength - 1], numFP[(int)cyclicPerm[wordLength - 1]], creal(fix(buffMat)), cimag(fix(buffMat))); 
		numFP[(int)cyclicPerm[wordLength - 1]]++;
		for (int j = 0; j < wordLength; j++){//Get the inverse of the current cyclic perm
			cyclicPerm[j] = (specialWord[(j + i) % wordLength] + 2) % 4;
//			printf("%d", cyclicPerm[j]);
		}
	//	printf("\n");
		makeWord(buffMat, gens, cyclicPerm, wordLength);
		fixRep[(int)cyclicPerm[wordLength - 1] * sizeArr + numFP[(int)cyclicPerm[wordLength - 1]]] = fix(buffMat);
		//printf("fix: %+lf %+lf\n", creal(fix(buffMat)), cimag(fix(buffMat)));
//		printf("fp[%d][%d] = %lf + %lf\n",(int)cyclicPerm[wordLength - 1], numFP[(int)cyclicPerm[wordLength - 1]], creal(fix(buffMat)), cimag(fix(buffMat))); 
		numFP[(int)cyclicPerm[wordLength - 1]]++;
	}

	//Now do the same for the abAB special word	
	char abABWord[4];
	abABWord[0] = 0;
	abABWord[1] = 1;
	abABWord[2] = 2;
	abABWord[3] = 3;
	wordLength = 4;
	cyclicPerm = (char*) calloc(wordLength, sizeof(char));
	for (int i = 0; i < wordLength; i++){
		for (int j = 0; j < wordLength; j++){
			cyclicPerm[j] = abABWord[(j + i) % wordLength];
//			printf("%d", cyclicPerm[j]);
		}
		makeWord(buffMat, gens, cyclicPerm, wordLength);
//		printf("\n");
		
		fixRep[(int)cyclicPerm[wordLength - 1] * sizeArr + numFP[(int)cyclicPerm[wordLength - 1]]] = fix(buffMat);
		//printf("fix: %+lf %+lf\n", creal(fix(buffMat)), cimag(fix(buffMat)));
//		printf("fp[%d][%d] = %lf + %lf\n",(int)cyclicPerm[wordLength - 1], numFP[(int)cyclicPerm[wordLength - 1]], creal(fix(buffMat)), cimag(fix(buffMat))); 
		numFP[(int)cyclicPerm[wordLength - 1]]++;
		for (int j = 0; j < wordLength; j++){//Get the inverse of the current cyclic perm
			cyclicPerm[j] = (abABWord[(j + i) % wordLength] + 2) % 4;
//			printf("%d", cyclicPerm[j]);
		}
//		printf("\n");
		makeWord(buffMat, gens, cyclicPerm, wordLength);
		fixRep[(int)cyclicPerm[wordLength - 1] * sizeArr + numFP[(int)cyclicPerm[wordLength - 1]]] = fix(buffMat);
		//printf("fix: %+lf %+lf\n", creal(fix(buffMat)), cimag(fix(buffMat)));
//		printf("fp[%d][%d] = %lf + %lf\n",(int)cyclicPerm[wordLength - 1], numFP[(int)cyclicPerm[wordLength - 1]], creal(fix(buffMat)), cimag(fix(buffMat))); 
		numFP[(int)cyclicPerm[wordLength - 1]]++;
	}

	/*
	for (int i = 0; i < 4; i++){
		for (int j = 0; j < numFP[i]; j++){
			printf("fp[%d][%d] : %lf + %lf\n", i, j, creal(fixRep[i * sizeArr + j]), cimag(fixRep[i * sizeArr + j])); 
		}
	}
	printf("%d\n", numFP[0]);
	printf("%d\n", numFP[1]);
	printf("%d\n", numFP[2]);
	printf("%d\n", numFP[3]);
	
	//bAba
	matmul(buff_gen_b, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][0] = fix(buff_out0);

	//aBa
	matmul(buff_gen_a, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][1] = fix(buff_out1);

	//Baa
	matmul(buff_gen_B, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	fixRep[0][2] = fix(buff_out1);

	//BAba
	matmul(buff_gen_B, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_a, buff_out0);
	fixRep[0][3] = fix(buff_out0);

	//ABab
	matmul(buff_gen_A, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][0] = fix(buff_out0);

	//AAb
	matmul(buff_gen_A, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	fixRep[1][1] = fix(buff_out1);

	//aBAb
	matmul(buff_gen_a, buff_gen_B, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_b, buff_out0);
	fixRep[1][2] = fix(buff_out0);

	//BabA
	matmul(buff_gen_B, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_b, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][0] = fix(buff_out0);

	//AbA
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][1] = fix(buff_out1);

	//bAA
	matmul(buff_gen_b, buff_gen_A, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	fixRep[2][2] = fix(buff_out1);

	//baBA
	matmul(buff_gen_b, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	matmul(buff_out1, buff_gen_A, buff_out0);
	fixRep[2][3] = fix(buff_out0);

	//abAB
	matmul(buff_gen_a, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_A, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][0] = fix(buff_out0);

	//aaB
	matmul(buff_gen_a, buff_gen_a, buff_out0);
	matmul(buff_out0, buff_gen_B, buff_out1);
	fixRep[3][1] = fix(buff_out1);

	//AbaB
	matmul(buff_gen_A, buff_gen_b, buff_out0);
	matmul(buff_out0, buff_gen_a, buff_out1);
	matmul(buff_out1, buff_gen_B, buff_out0);
	fixRep[3][2] = fix(buff_out0);
	*/

}

void computeCycles(double complex begpt[4], double complex endpt[4], double complex* gens){

	double complex buff_gen0[2][2];
	double complex buff_gen1[2][2];
	double complex buff_gen2[2][2];
	double complex buff_gen3[2][2];
	double complex buff_out0[2][2];
	double complex buff_out1[2][2];
	//Copy gens to buffers (since I couldn't find a clean way to matmul)
	matrix3dto2D(gens, buff_gen0, 0);
	matrix3dto2D(gens, buff_gen1, 1);
	matrix3dto2D(gens, buff_gen2, 2);
	matrix3dto2D(gens, buff_gen3, 3);

	//Compute the fix points of all right most turns 
	matmul(buff_gen3, buff_gen2, buff_out0); 
	matmul(buff_out0, buff_gen1, buff_out1); 
	matmul(buff_out1, buff_gen0, buff_out0); 
	endpt[0] = fix(buff_out0);

	matmul(buff_gen0, buff_gen3, buff_out0); 
	matmul(buff_out0, buff_gen2, buff_out1); 
	matmul(buff_out1, buff_gen1, buff_out0); 
	endpt[1] = fix(buff_out0);

	matmul(buff_gen1, buff_gen0, buff_out0); 
	matmul(buff_out0, buff_gen3, buff_out1); 
	matmul(buff_out1, buff_gen2, buff_out0); 
	endpt[2] = fix(buff_out0);

	matmul(buff_gen2, buff_gen1, buff_out0); 
	matmul(buff_out0, buff_gen0, buff_out1); 
	matmul(buff_out1, buff_gen3, buff_out0); 
	endpt[3] = fix(buff_out0);

	//Compute the fix points of all right most turns 
	matmul(buff_gen1, buff_gen2, buff_out0); 
	matmul(buff_out0, buff_gen3, buff_out1); 
	matmul(buff_out1, buff_gen0, buff_out0); 
	begpt[0] = fix(buff_out0);

	matmul(buff_gen2, buff_gen3, buff_out0); 
	matmul(buff_out0, buff_gen0, buff_out1); 
	matmul(buff_out1, buff_gen1, buff_out0); 
	begpt[1] = fix(buff_out0);

	matmul(buff_gen3, buff_gen0, buff_out0); 
	matmul(buff_out0, buff_gen1, buff_out1); 
	matmul(buff_out1, buff_gen2, buff_out0); 
	begpt[2] = fix(buff_out0);

	matmul(buff_gen0, buff_gen1, buff_out0); 
	matmul(buff_out0, buff_gen2, buff_out1); 
	matmul(buff_out1, buff_gen3, buff_out0); 
	begpt[3] = fix(buff_out0);

}

