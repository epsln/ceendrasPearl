#ifndef treeExploration
#define treeExploration

void goForward(int *lev, int* tag, double complex word[1000][2][2], double complex gens[4][2][2]);
void goBackwards(int *lev);
int availableTurn(int *lev, int* tag);
void turnForward(int *lev, int tag[1000], double complex word[1000][2][2], double complex gens[4][2][2]);
void computeDepthFirst(double* PARAMS, double complex ta, double complex tb, float*** imgArr, int numIm);
int branchTerm(double* PARAMS, double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex fixRep[4][3], double complex word[1000][2][2], float*** imgArr);
int branchTermEpsi(double* PARAMS, double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex fixRep[4][3], double complex word[1000][2][2], float*** imgArr);
#endif
