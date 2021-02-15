#ifndef treeExploration
#define treeExploration

void goForward(int *lev, int* tag, int* state, int FSA[19][4], mpc_t* word, mpc_t* gens);
void goBackwards(int *lev);
int availableTurn(int *lev, int *tag, int* state, int FSA[19][4]);
void turnForward(int *lev, int *tag,int* state, int FSA[19][4], mpc_t* word, mpc_t* gens);

int branchTermRepetends(int lev, int* tag, mpc_t fixRep[4][4], mpc_t* word, image_t* img);

void computeDepthFirst(double complex* gens, image_t* img, int numIm);
#endif
