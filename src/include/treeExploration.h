#ifndef treeExploration
#define treeExploration

void goForward(int *lev, int* tag, double complex word[1000][2][2], double complex gens[4][2][2]);
void goBackwards(int *lev);
int availableTurn(int *lev, int* tag);
void turnForward(int *lev, int tag[1000], double complex word[1000][2][2], double complex gens[4][2][2]);

int branchTermEpsi(double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex word[1000][2][2], image_t* img);//TODO: Remove this (deprecated)
int branchTermRepetends(double complex* oldPoint, int lev, int* tag, double complex fixRep[4][4], double complex word[1000][2][2], image_t* img);

void computeDepthFirst(double complex gens[4][2][2], image_t* img, int numIm);
#endif
