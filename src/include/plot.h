#ifndef plot
#define plot

typedef struct image{
	//Image characteristics
	int w;
	int h;
	int line;
	int antialiasingPow;
	int bitwise; //Flag to use the bitwise ops
	int debug;


	//Computation parameters
	int levmax;
	double bounds;
	double epsi;

	char* filename;

	int *pointArr;//1D array to store the endpoints of the exploration

	//Experimental !	
	//Use an int array to store point flags
	//One int_64 (long int) holds flags for a line of 64 points 
	//We then use bitwise ops to do a fast downscale
	//...At least that's the theory
	unsigned long long int *bitArray;//1D array to store the endpoints of the exploration
	

}image_t;


int checkBoundaries(int x, int y, image_t* img);
void plotLineLow(int x0,int y0, int x1,int y1, image_t* img);
void plotLineHigh(int x0,int y0, int x1,int y1, image_t* img );
void point(int x,int y, image_t* img);
void line(int x0,int y0, int x1,int y1,  image_t* img );
void antialiasing(image_t* img, unsigned char* output);
void saveArrayAsBMP(image_t* img);

#endif
