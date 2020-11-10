#ifndef plot
#define plot

typedef struct image{
	int w;
	int h;
	int levmax;
	int line;
	int antialiasingPow;
	int debug;

	double bounds;
	double epsi;

	char* filename;

	int *pointArr;//2D array to store the endpoints of the exploration
	

}image_t;


int checkBoundaries(int x, int y, image_t* img);
void plotLineLow(int x0,int y0, int x1,int y1, image_t* img);
void plotLineHigh(int x0,int y0, int x1,int y1, image_t* img );
void point(int x,int y, image_t* img);
void line(int x0,int y0, int x1,int y1,  image_t* img );
void antialiasing(image_t* img, float*** output);
void saveArrayAsBMP(image_t* img);

#endif
