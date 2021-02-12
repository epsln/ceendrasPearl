#include <stdio.h>
#include <time.h>
#include <math.h>

#include "include/progressBar.h"

void pBarAnim(int numImg, int totalImg, double timeArray[totalImg]){
	//One pbar that shows everything that needs to be done in order to finish the program
	//Shows time per image
	//Shows ETA
	//It's not a progress bar atm, but will be sometime

	double avgTimeImage = 0;
	double diff = 0;
	//Time Array is an overkill method to estimate time per image
	//Each time an image is done, the timestamp is saved at the numImg position in the array
	//Then we can do some averaging to get the mean time per image and use this for our estimate
	time_t rawtime;
	struct tm * eta;

	time ( &rawtime );

	timeArray[numImg] = (double)clock()/CLOCKS_PER_SEC;
	avgTimeImage = timeArray[0];

	for (int i = 1; i <= numImg; i++){
		avgTimeImage += timeArray[i] - timeArray[i - 1];	
	}
	avgTimeImage /= numImg + 1; 
	rawtime = rawtime + avgTimeImage * (totalImg - numImg);
	eta = localtime ( &rawtime );

	diff = timeArray[numImg] - timeArray[numImg - 1];

	printf("Compute time for last image: %d:%d:%d:%d\n", (int) diff / 3600, (int) diff /  60, (int) diff % 60, (int) trunc(diff * 1000));
	printf("Average compute time: %d:%d:%d:%d\n", (int) avgTimeImage / 3600, (int) avgTimeImage / 60 , (int) avgTimeImage % 60, (int) trunc(avgTimeImage * 1000));
	printf("ETA: %d/%d/%d %d:%d:%d\n",eta->tm_mday, eta->tm_mon + 1, eta->tm_year + 1900, eta->tm_hour, eta->tm_min, eta->tm_sec);
	//A quick hacky progess bar 
	printf("[");
	for (int i = 0; i < (float)numImg/(totalImg) * 50; i++){
		printf("=");
	}
	printf(">");
	for (int i = (float)numImg/(totalImg) * 50; i < 50; i++){
		printf(" ");
	}
	printf("] %.2f \% \n\n", (float)numImg/totalImg * 100 );
}
