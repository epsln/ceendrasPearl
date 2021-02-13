#include <stdio.h>
#include <time.h>
#include <math.h>

#include "include/progressBar.h"

void pBarAnim(int numImg, int totalImg, double timeArray[10]){
	//One pbar that shows everything that needs to be done in order to finish the program
	//Shows time per image
	//Shows ETA
	//It's not a progress bar atm, but will be sometime

	double avgTimeImage = 0;
	double diff = 0;

	//Each time an image is done, the number of sec since prog launch is saved at the numImg position in the array
	//Then we can do some floating average to get the mean time per image and use this for our estimate
	time_t rawtime;
	struct tm * eta;

	//Ge the current time
	time ( &rawtime );


	//Move all values one idx back
	if(numImg > 10){
		for (int i = 0; i < fmin(10, numImg); i++){
			timeArray[i] = timeArray[i + 1]; 
		}
	}

	//Add at current time at the end
	//Divide clock by clock per sec to get time in sec
	timeArray[(int)fmin(9, numImg)] = (double)clock()/CLOCKS_PER_SEC; 

	//Rolling average 
	for (int i = 1; i < fmin(10, numImg); i++){
		avgTimeImage += timeArray[i] - timeArray[i - 1];	
	}
	avgTimeImage /= fmin(numImg + 1, 10); 

	//Add the average time * num of image left to get ETA
	rawtime = rawtime + avgTimeImage * (totalImg - numImg);

	//Get that ETA in a nice format 
	eta = localtime ( &rawtime );

	diff = timeArray[(int)fmin(9, numImg)] - timeArray[(int)fmin(8, numImg - 1)];

	printf("Compute time for last image: %d:%d:%d:%d\n", (int) diff / 3600, (int) diff /  60, (int) diff % 60, (int) trunc(diff * 1000));
	printf("Average compute time: %d:%d:%d:%d\n", (int) avgTimeImage / 3600, (int) avgTimeImage / 60 , (int) avgTimeImage % 60, (int) trunc(avgTimeImage * 1000));
	if (numImg == totalImg - 1){
		diff = timeArray[numImg] - timeArray[0];
		printf("Total Compute Time: %d:%d:%d:%d\n", (int) diff / 3600, (int) diff /  60, (int) diff % 60, (int) trunc(diff * 1000));
		printf("Thanks for your patience ! Love you bye\n");
	}
	else	
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
