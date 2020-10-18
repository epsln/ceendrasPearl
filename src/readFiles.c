#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>


void readPoints(double *pointsList[2]){
	FILE* pointsX = NULL;
	FILE* pointsY = NULL;

	pointsX = fopen("../jungPtsX.txt", "r");
	pointsY = fopen("../jungPtsY.txt", "r");

	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	double x;
	int i = 0;
	if(pointsX == NULL || pointsY == NULL){
		printf("Couldn't open documents !\n");
		exit(-1);
	}

	while ((read = getline(&line, &len, pointsX)) != -1) {
		x = atof(line);
		pointsList[i][0] = x;	
		i+=1;
	}
	i = 0;
	while ((read = getline(&line, &len, pointsY)) != -1) {
		x = atof(line);
		pointsList[i][1] = x;	
		i+=1;
	}
	fclose(pointsX);
	fclose(pointsY);
}

void readConf(int *RESX, int *RESY, int *NPOINTS, int *JITTER, double *JITTER_B, int *GRAY, int *OPENCL_ITER, int *MAXITER,int *MAXITER_G, int *MAXITER_B, double *RED_C, double *GRE_C, double *BLU_C, char filename[256], int *RAND, char kernelFilename[256]){
	FILE* conf = NULL;

	conf = fopen("./config.cfg", "r");
	char * line = NULL;
	char buff[256] = "";
	int count = 0;

	size_t len = 0;
	ssize_t read;

	if (conf == NULL){
		printf("Could not open config file !\n");
		exit(-1);
	}

	while ((read = getline(&line, &len, conf)) != -1) {

		switch(count){
			case 1:
				*RAND = atoi(line);
				break;
			case 3:
				*NPOINTS = atoi(line);
				break;
			case 5:
				*RESX = atoi(line);
				break;
			case 7:
				*RESY = atoi(line);
				break;	
			case 9:
				*JITTER = atoi(line);
				break;	
			case 11:
				*JITTER_B = atof(line);
				break;	
			case 13:
				*GRAY = atoi(line);
				break;	
			case 15:
				*OPENCL_ITER = atoi(line);
				break;	
			case 17:
				*MAXITER = atoi(line);
				break;	
			case 19:
				*MAXITER_G = atoi(line);
				break;	
			case 21:
				*MAXITER_B = atoi(line);
				break;	
			case 23:
				*RED_C = atof(line);
				break;	
			case 25:
				*GRE_C = atof(line);
				break;	
			case 27:
				*BLU_C = atof(line);
				break;	
			case 29:
				strcpy(buff, line);
				//Remove the tab that is on first pos
				for (int i = 0; i < 255; i++){
					filename[i] = buff[i+1];

				}
				strtok(filename, "\n");//Remove newline
				break;
			case 31:
				strcpy(buff, line);
				//Remove the tab that is on first pos
				for (int i = 0; i < 255; i++){
					kernelFilename[i] = buff[i+1];

				}
				strtok(kernelFilename, "\n");//Remove newline
				break;



		}
		count++;

	}


	fclose(conf);	
}
