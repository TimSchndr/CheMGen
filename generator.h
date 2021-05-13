#define matrix(i, j) (adjacencyMatrix[i * numberOfAtoms + j])
#define MIN_VALUE(a, b) ((a) < (b) ? a : b)
#define TRUE 1
#define FALSE 0

void setAtomValency();
int getAtomValency(char atom);
char* defineAtomName();
int checkAtoms();
void makeBond(int i, int j, int bondSum);
void saveGraph();
void writeData(char* input);
void dfs(int start);
int checkConnectivity();
int getMinimumValue(int a, int b, int c);
int freeAll();

int* getComponents(char* formula);
int* defineAlphabet(char* string);
int calcDoubleBounds(const int* components);
int* addElements(const int* components);