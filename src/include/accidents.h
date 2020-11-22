#ifndef accidents
#define accidents 

//This is all function to find "accidents": generator who have a parabolic combinations 
//For example, tr(aB) = tr(b) = 2
//Here aB can be designed as the irreducible fraction 1/2
//Basically, we implements everything and anything in the Ch. 9 :)

double complex traceRecursion(int p, int q, double complex ta, double complex tB, double complex taB);
double complex findMu(int denom);//Slight modification, here we don't want to plot the boundary, but get the value of mu to plot the associated "accident"
double newtPQ(int p, int q, int denom);//Gives out the next item in the farray sequence with denum being the maximum denominator
double newtonSolver(double complex initGuess);//Newton solve our way into a mu existence


#endif
