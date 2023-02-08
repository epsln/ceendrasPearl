#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <complex.h>

#include "include/complexMath.h"
#include "include/readFiles.h"
#include "include/plot.h"

void readPoints(double complex pointsList[256][2], image_t *img){
	FILE* fp = NULL;

	fp = fopen("./pointList.cfg", "r");

	char line[80];

	double a, b, c, d;
	int i = 0;
	if(fp == NULL){
		printf("Couldn't open points file!\n");
		exit(-1);
	}
	while (fgets(line, 80, fp) != NULL) {
		if (sscanf(line, "%lf %lfi, %lf %lfi\n", &a, &b, &c, &d) == 4){
			pointsList[i][0] = a + I * b;	
			pointsList[i][1] = c + I * d;	
		}
		if (strcmp(line, "RANDOM") == 1){
			pointsList[i][0] = randomComplex(-img -> bounds, img -> bounds);
			pointsList[i][1] = randomComplex(-img -> bounds, img -> bounds);
		}
		if (strcmp(line, "RANDOM FIX DIST") == 0 && i == 0){
			pointsList[i][0] = randomComplex(-img -> bounds, img -> bounds);
			pointsList[i][1] = randomComplex(-img -> bounds, img -> bounds);
		}
		if (strcmp(line, "RANDOM FIX DIST") == 0 && i != 0){
			pointsList[i][0] = randomComplexFixDist(pointsList[i-1][0], a, img -> bounds);
			pointsList[i][1] = randomComplexFixDist(pointsList[i-1][1], a, img -> bounds);
		}
		//if (sscanf(line, "%d/%d", &a, &b) == 2)
		//		getTraceFromFract(pointsList[i], (ratio) {a, b});	
		i+=1;
	}
	fclose(fp);
}

/*TODO: Finish up, and modify every important thing from config file
  void readConf(char filename[256], image_t *img){
  FILE* fp = NULL;

  fp = fopen("./config.cfg", "r");
  char* line = NULL;
  char buff[256] = "";
  double a, b;

  if (fp == NULL){
  printf("Could not open config file !\n");
  exit(-1);
  }

  while (fgets(line, 80, fp) != NULL) {
  if (sscanf(line, "%lf + %lfi", &a, &b) == 2)
  pointsList[i] = a + I * b;	
  if (strcmp(line, "RANDOM") == 0)
  pointsList[i] = randomComplex(-img -> bounds, img -> bounds);
  if (strcmp(line, "RANDOM FIX DIST") == 0 && i == 0)
  pointsList[i] = randomComplex(-img -> bounds, img -> bounds);
  if (strcmp(line, "RANDOM FIX DIST") == 0 && i != 0)
  pointsList[i] = randomComplexFixDist(pointsList[i-1], a, img -> bounds);
  i+=1;
  }
  fclose(fp);
  }
  */
