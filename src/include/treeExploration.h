#ifndef treeExploration
#define treeExploration

void goForward(int *lev, int* tag, double complex word[1000][2][2], double complex gens[4][2][2]);
void goBackwards(int *lev);
int availableTurn(int *lev, int* tag);
void turnForward(int *lev, int tag[1000], double complex word[1000][2][2], double complex gens[4][2][2]);
void computeDepthFirst(double complex ta, double complex tb,  double complex tab, image_t* img, int numIm);
int branchTermEpsi(double complex* oldPoint, int lev, int* tag, double complex endpt[4], double complex word[1000][2][2], image_t* img);
int branchTermRepetends(double complex* oldPoint, int lev, int* tag, double complex fixRep[4][4], double complex word[1000][2][2], image_t* img);
#endif
