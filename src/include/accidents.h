#ifndef accidents
#define accidents 

//This is all function to find "accidents": generator who have a parabolic combinations 
//For example, tr(aB) = tr(b) = 2
//Here aB can be designed as the irreducible fraction 1/2
//Basically, we implements everything and anything in the Ch. 9 :)
//It would seem that the nextPQ algo given by the book doesnt work, so we'll work around...
typedef struct {
	int p, q;
}ratio;

double complex tracePoly(ratio fraction, double complex ta, double complex tB, double complex taB);
double complex traceEqn(ratio fraction, double complex mu);
void nextPQ(int *p, int *q, int denom);//Gives out the next item in the farray sequence with denum being the maximum denominator
void makeFareySeq(int denom, ratio* fareyArr);//Gives out the next item in the farray sequence with denum being the maximum denominator
void newtonSolver(double complex *pZ0, ratio fraction);//Newton solve our way into a mu existence


#endif
