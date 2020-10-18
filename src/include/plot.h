#ifndef plot
#define plot
int checkBoundaries(int x, int y, int maxW, int maxH);
void plotLineLow(int x0,int y0, int x1,int y1, float*** imgArr);
void plotLineHigh(int x0,int y0, int x1,int y1, float*** imgArr);
void line(int x0,int y0, int x1,int y1, float*** imgArr, int LINE, int WIDTH, int HEIGHT);
void saveArrayAsBMP(float*** imgArr, char* filename, double* PARAMS);

#endif
