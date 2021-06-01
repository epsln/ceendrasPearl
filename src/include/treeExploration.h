#ifndef treeExploration
#define treeExploration

typedef struct{
	double complex* gens;
	image_t* img;
	int numIm;
	int numBranch;
	char* specialWord;
	double complex* fixRep;
	int wordLength;
	int numFP[4];
}dfsArgs;

void goForward(int *lev, int* tag, int* state, int FSA[19][4], double complex* word, double complex* gens);
void goBackwards(int *lev);
int availableTurn(int *lev, int *tag, int* state, int FSA[19][4]);
void turnForward(int *lev, int *tag,int* state, int FSA[19][4], double complex* word, double complex* gens);

int branchTermRepetends(int lev, int* tag, double complex* fixRep, int wordLength, int numFP[4], double complex* word, image_t* img);

//void computeDepthFirst(double complex* gens, image_t* img, int numIm);
void *computeDepthFirst(void *_args);
#endif
