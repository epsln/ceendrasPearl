#ifndef treeExploration
#define treeExploration

void goForward(int *lev, int* tag, int* state, int FSA[19][4], double complex* word, double complex* gens);
void goBackwards(int *lev);
int availableTurn(int *lev, int *tag, int* state, int FSA[19][4]);
void turnForward(int *lev, int *tag,int* state, int FSA[19][4], double complex* word, double complex* gens);

int branchTermRepetends(int lev, int* tag, double complex fixRep[4][4], double complex* word, image_t* img);

void computeDepthFirst(double complex* gens, image_t* img, int numIm);

void goForwardQuat(int *lev, int* tag, int* state, int FSA[19][4], quat_t* word, quat_t* gens);
void turnForwardQuat(int *lev, int *tag,int* state, int FSA[19][4], quat_t* word, quat_t* gens);
int branchTermRepetendsQuat(int lev, int* tag, quat_t fixRep[4][4], quat_t* word, image_t* img);
void computeDepthFirstQuat(quat_t* gens, image_t* img, int numIm);
#endif
