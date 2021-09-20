#define matrix(i, j) (adjacencyMatrix[i * numberOfnonHAtoms + j])
#define MIN_VALUE(a, b) ((a) < (b) ? a : b)
#define TRUE 1
#define FALSE 0

void setAtomValence();
int getAtomValence(char atom);
void defineAtomName(int* components, int differentAtoms);
int checkAtoms(int** partitionList, int differentAtoms);
int checkAtoms_1(int differentAtoms);
void makeBond(int i, int j, int bondSum, int** partition, int differentAtoms);
void makeBond_1(int i, int j, int bondSum, int differentAtoms);
void dfs(int start, int* visitedNodes);
int checkConnectivity();
int shortestCycle();
int getNextComposition(int numberOfH, int length);
int checkComposition(int differentAtoms);
int* allocatePartition(int length);
void setPartitionToZero(int length, int* partition);
int nextPartition(int number, int* array);
int nextAtomPartition(int numberOfH, int arrayLength, int* partition, char atomType);
int checkAtomPartition(int number, int* partition);
int getPartitionCount(int index, int** partition);
void writeSDFformat();
int freeAll();
long* calculateMorganNumbers();
int morganNumbersExist(long* molecule);
void saveMorganNumbers(long* molecule);
int countDigits(long* array, int length);

int* getComponents(char* formula);
int* defineAlphabet(char* string);
int calcDoubleBonds(const int* components);
int countAtomTypes(int* components);
int* addElements(const int* components);

void quickSort(long* left, long* right);