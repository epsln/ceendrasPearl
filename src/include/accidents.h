#ifndef accidents
#define accidents 

//This is all function to find "accidents": generator who have a parabolic combinations 
//For example, tr(aB) = tr(b) = 2
//Here aB can be designed as the irreducible fraction 1/2
//Basically, we implements everything and anything in the Ch. 9 :)
typedef struct {
	long long int p, q;
}ratio;

ratio simplify_fract(ratio a);

double complex tracePoly(ratio fraction, double complex ta, double complex tB, double complex taB);
double complex traceEqn(ratio fraction, double complex mu);
void makeFareySeq(int denom, ratio* fareyArr);//Populates a ratio array with all farey sequence with maximum denominator denum
void makeFiboSeq(int lengthAnim, ratio* fareyArr);//Populates a ratio array with the Fibonacci Sequence ratio (F_{n+1}/F_n)
void makePiSeq(int lengthAnim, ratio* fareyArr);//Populates a ratio array with the Fibonacci Sequence ratio (F_{n+1}/F_n)
int makeContinuedFraction(int lengthAnim, double real, ratio* fareyArr);//Populates a ratio array with the Fibonacci Sequence ratio (F_{n+1}/F_n)

void getTraceFromFract(double complex *pz0, ratio fraction);
void getSpecialWordFromFract(ratio fraction, char* specialWord);//Get the special word out of the fraction representation

void newtonSolver(double complex *pZ0, ratio fraction);//Newton solve our way into a mu existence


#endif
