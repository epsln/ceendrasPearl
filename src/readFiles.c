#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#include "include/plot.h"

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

float getNumFromLine(char* line){
	//Helper function to get a number from a line in a config file
	//The number must be the last thing on a line
	//STRTOK BITCH
	int idx = 0;
	char numArr[256]; 
	for (int i = 0; i < 256; i++){
		if (line[i] == " ")
			idx = i;
	}
 	while (line[idx] != '\n' || line[idx] != '\0'){
		numArr[idx] = line[idx];
		idx++;
	}	
	numArr[idx] = '\0';
	return atoi(numArr);
}

void readConf(image_t *pImg){
	FILE* conf = NULL;

	conf = fopen("./params.cfg", "r");
	char * line = NULL;
	char buff[256] = "";
	int linecount = 0;

	size_t len = 0;
	ssize_t read;

	if (conf == NULL){
		printf("Could not open parameter file ! Is params.cfg in same folder as the executable ?\n");
		exit(-2);
	}

	while ((read = getline(&line, &len, conf)) != -1) {
		if (line[0] = "#")
			linecount++;
			continue;
		if(strcmp(strtok(line, " "), "ANTIALPOW"))
				pImg -> antialiasingPow =atoi(strtok(NULL, " ")) ;
		else if (strcmp(strtok(line, " "), "WIDTH"))
				pImg -> w = atoi(strtok(NULL, " "));
		else if (strcmp(strtok(line, " "), "HEIGHT"))
				pImg -> h = atoi(strtok(NULL, " "));
		else if (strcmp(strtok(line, " "), "BOUNDS"))
				pImg -> bounds =atoi(strtok(NULL, " "));
		else if (strcmp(strtok(line, " "), "EPSI"))
				pImg -> epsi =atoi(strtok(NULL, " ")); 
		else if (strcmp(strtok(line, " "), "MAXLVL"))
				pImg -> maxword =atoi(strtok(NULL, " "));
		else if (strcmp(strtok(line, " "), "BITWISE"))
				pImg -> bitwise =atoi(strtok(NULL, " "));
	fclose(conf);	
}
