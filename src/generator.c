#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>

#include "generator.h"

// X stands for Chlorine
// B stands for Bromine, *not* Boron
#define ORGANIC_SUBSET "BbCcFfHhIiNnOoPpSsXx"
#define ELEMENTS "BCFHINOPSX"


//general information
int atomValency[26] = { 0 };

int* components;
int numberOfAtoms;
int* pAtoms = &numberOfAtoms;
int numberOfBonds;
int* pBonds = &numberOfBonds;
char* atomName;

int counter = 0;
int* pCounter = &counter;

//the graph properties
//this is a symmetric Laplacian matrix
int* adjacencyMatrix;
int* visitedNodes;

FILE* outputFile;

/*
* Driver Code
*/
void main(int argc, char* argv[]) {

	//if no file is specified, write to stdout
	if (argc == 2) {
		outputFile = stdout;
	}
	else if (argc == 3) {
		outputFile = fopen(argv[2], "w");
	}
	else {
		fprintf(stderr, "Usage: generator <forumla> <outputFile>");
	}

	if (outputFile) {

		//get the components of the formula
		components = getComponents(argv[1]);

		//get the number of bonds
		numberOfBonds = *pAtoms - 1 + calcDoubleBounds(components);

		//define list of atoms
		atomName = defineAtomName();

		printf("number of atoms (non-H): %d\nnumber of bonds: %d\n\n", *pAtoms, *pBonds);

		//write general data to file
		writeData(argv[1]);

		//allocation of adjacencyMatrix
		adjacencyMatrix = calloc(*pAtoms * *pAtoms, sizeof(adjacencyMatrix));

		if (adjacencyMatrix) {

			//set valency of atoms
			setAtomValency();
			//start filling matrix
			makeBond(0, 1, 0);

		}
		else {
			fprintf(stderr, "Allocation failed for adjacencyMatrix in main().");
			freeAll();
			exit(0);
		}

		//only if there was an output file close it
		if (outputFile != stdout) {
			fclose(outputFile);
		}

		printf("number of valid graphs: %d\n", *pCounter);
	}
	else {
		fprintf(stderr, "Could not open the specified file.");
		freeAll();
		exit(0);
	}
}

/*
* A method that sets the number of valencies for each atom.
*
* B - Brome
* X - Chlorine
*/
void setAtomValency() {

	// (int) 'A' == 65
	atomValency['B' - 65] = 1;
	atomValency['C' - 65] = 4;
	atomValency['F' - 65] = 1;
	atomValency['H' - 65] = 1;
	atomValency['I' - 65] = 1;
	atomValency['N' - 65] = 3;
	atomValency['O' - 65] = 2;
	atomValency['P' - 65] = 2;
	atomValency['S' - 65] = 2;
	atomValency['X' - 65] = 1;

}

/*
* A method that returns the valency of the given atom
*
* char atom - the atom abbreviation for which the valency is required
* return - the atoms' valency
*/
int getAtomValency(char atom) {

	return atomValency[atom - 65];

}

/*
* A method that creates an array with the atom identifier for each atom in the formula
*
* return - the pointer to the array with the atom names
*/

char* defineAtomName() {

	char* atomList = calloc(*pAtoms + 1, sizeof(atomList));

	if (atomList) {
		int index = 0;

		//for all components not 'H' fill the array with the atom names
		for (unsigned int i = 0; i < 26; i++) {
			for (int j = 0; j < components[i]; j++) {

				if (i + 65 != 'H') {
					atomList[index] = i + 65;
					index++;
				}
			}
		}
	}
	else {
		fprintf(stderr, "Allocation failed in defineAtomName().");
		freeAll();
		exit(0);
	}

	return atomList;
}

/*
* A Method that checks if the valency of all atoms in the matrix is correct in comparrison
* to their natural maximum valency.
*
* return - an int telling if all atoms are correct according to valency
*/
int checkAtoms() {

	for (int i = 0; i < *pAtoms; i++) {
		if (matrix(i, i) > getAtomValency(atomName[i])) {
			return FALSE;
		}
	}

	return TRUE;
}

/*
* A method that fills the adjacencyMatrix with bonds recursively.
*
* int i - the row index
* int j - the column index
* int bondSum - the number of bonds already spent
*/
void makeBond(int i, int j, int bondSum) {

	if (*pBonds == bondSum) {
		// all bonds have been placed

		if (checkAtoms()) {
			// check if all atoms are ok

			if (checkConnectivity()) {
				//check if the graph is connected

				//save the graph
				saveGraph();

				counter++;
				return;
			}
			return;
		}
		return;
	}

	//check if the end of the matrix is reached
	if ((i * *pAtoms + j) >= (*pAtoms * *pAtoms)) {
		return;
	}

	int rememberValue = matrix(i, j);

	//add a new bond between i and j
	for (int b = 0; b <= MIN_VALUE(getAtomValency(atomName[i]), getAtomValency(atomName[j])); b++) {

		matrix(i, j) = b;
		matrix(i, i) += b;
		matrix(j, j) += b;

		if (j < *pAtoms - 1) {
			makeBond(i, j + 1, bondSum + b);
		}
		else {
			makeBond(i + 1, i + 2, bondSum + b);
		}
		matrix(i, i) -= b;
		matrix(j, j) -= b;
	}

	matrix(i, j) = rememberValue;

}

/*
* A method that appends a graph to the file.
* The matrix is written as upper triangle matrix without the diagonal.
*/

void saveGraph() {

	for (int i = 0; i < *pAtoms; i++) {
		for (int j = i + 1; j < *pAtoms; j++) {
			fprintf(outputFile, "%d ", matrix(i, j));
		}
	}
	fprintf(outputFile, "\n");

}

/*
* A method that creates an new file in the given directory, initializing it with the
* neccessary data, such as molecular formula, number of Atoms, number of bonds.
*
* char* input - the formula read from cmd
*/

void writeData(char* input) {

	// write the formula to file
	fprintf(outputFile, "%s\n", input);
	// write the number of non-H atoms to file
	fprintf(outputFile, "%d\n", *pAtoms);

	//write each atom identifier in the correct order (order of internal ussage)
	for (unsigned int i = 0; i < strlen(atomName); i++) {
		fprintf(outputFile, "%c ", atomName[i]);
	}

	fprintf(outputFile, "\n");

}

/*
* A method that does a DFS on an adjacency matrix for given start node.
*
* int start - the node to start from
*/
void dfs(int start) {

	visitedNodes[start] = 1;

	//only if the start node has degree > 0 the Graph can be connected
	if (matrix(start, start)) {

		//search for the next adjacent node to start from (not on main diagonal)
		for (int i = 0; i < *pAtoms; i++) {
			//if vertex i has not been visited and there is an edge from start to i, do dfs
			if (!visitedNodes[i] && matrix(start, i)) {
				dfs(i);
			}
		}
	}
}

/*
* A method that checks if the underlying graph is connected.
*/

int checkConnectivity() {

	visitedNodes = calloc(*pAtoms, sizeof(visitedNodes));

	if (visitedNodes) {

		//start the DFS with node 0
		dfs(0);

		//check if all nodes have been visited
		for (int i = 0; i < *pAtoms; i++) {
			if (!visitedNodes[i]) {
				free(visitedNodes);
				return FALSE;
			}
		}

		free(visitedNodes);
		return TRUE;
	}
	else {
		printf("Allocation of 'visitedNodes' failed.");
		freeAll();
		exit(0);
	}
}

/*
* A method that calculates the minimum of three values.
*
* int a, b, c - the three values
* return - the minimum of three values
*/
int getMinimumValue(int a, int b, int c) {

	int temp;

	if (a < b) {
		temp = a;
	}
	else {
		temp = b;
	}

	if (c < temp) {
		return c;
	}
	else {
		return temp;
	}
}

/*
* Free all global variables.
* This method is called when an error occured during run.
*/
int freeAll() {

	if (components) {
		free(components);
	}
	if (atomName) {
		free(atomName);
	}
	if (adjacencyMatrix) {
		free(adjacencyMatrix);
	}
	if (visitedNodes) {
		free(visitedNodes);
	}

	return 1;
}
/************************************************************************************************************/
/******************************************* Begin of Formula Reader ****************************************/
/************************************************************************************************************/
/*
* A method that reads in a String from console (molecular formula) and
* counts the number of occurences of each element.
*
* return - an array containing the count for each letter in the input
*/

// TODO case abfangen, dass Nummer an Pos 0 oder Char nicht in Alphabet

int* getComponents(char* formula) {

	int amount[256] = { 0 };
	int* alphabet = defineAlphabet(ORGANIC_SUBSET);

	//check if allocation went well
	if (alphabet && amount) {

		int checkForElement = TRUE;
		int identifier = -1;
		unsigned int i = 0;

		//test if there is a digit at position 0
		if (isdigit(formula[0])) {
			fprintf(stderr, "There must be an atom identifier at the beginning.\n");
			freeAll();
			exit(0);
		}

		while (i < strlen(formula)) {

			// expect an atom identifier at the next position
			if (checkForElement) {

				//check if the read character is in the alphabet
				if (!alphabet[formula[i]]) {
					fprintf(stderr, "Character at pos %d is not part of the alphabet.\n", i);
					freeAll();
					exit(0);
				}

				identifier = formula[i];

				// next string element is a number, memorize the atom
				if (isdigit(formula[i + 1])) {
					checkForElement = FALSE;
					i++;
				}
				//next string element is a character or end of string
				else if (isalpha(formula[i + 1]) || formula[i + 1] == '\0') {
					amount[identifier]++;
					i++;
				}

			}
			// expect a number at the next position
			else {
				int j = i;
				int len = 0;
				//search for the end of the number
				while (isdigit(formula[j])) {
					j++;
					len++;
				}

				char* num = calloc(len + 1, sizeof * num);

				if (num && len > 0) {

					memcpy(num, &formula[i], len + 1);

					//the number behind the current atom
					int absolut = atoi(num);
					amount[identifier] += absolut;
					free(num);
				}
				else {
					fprintf(stderr, "Allocation of num in getComponents() failed\n");
					freeAll();
					exit(0);
				}
				//jump behind the read number in the input
				i += len;
				checkForElement = TRUE;
			}
		}
	}
	free(alphabet);



	//add upper_case and lower_case counts to one variable
	int* elementCounts = addElements(amount);

	return elementCounts;
}

/*
* A method that computes the number of Double Bond Equivalents on basis of the given input.
*
* int* elementCounts - the count for each upper_case character in ASCII (length 26)
*
* return - the number of DBE for given input
*/

int calcDoubleBounds(const int* elementCounts) {


	//total number of halogens in the input
	int halogenSum = 0;
	char halogens[5] = "BFIX";

	for (unsigned int i = 0; i < strlen(halogens); i++) {
		halogenSum += elementCounts[halogens[i] - 65];
	}

	//the number of H, N and halogens combined according to DBE formula
	double negativeTerm = 0;
	negativeTerm = 0.5 * ((double)elementCounts['H' - 65] - (double)elementCounts['N' - 65] + (double)halogenSum);
	int dbe = 0;

	//check if negativeTerm is an int
	if (negativeTerm == (int)negativeTerm) {

		dbe = elementCounts['C' - 65] + 1 - negativeTerm;
	}
	else {
		fprintf(stderr, "The input was not correct with respect to DBE calculation.");
		freeAll();
		exit(0);
	}
	return dbe;
}

/*
* A method to define an ASCII based alphabet on basis of a given string.
*
* char* string - the string to use as reference
*
* return - an array filled with 1 for each character in the string
*/

int* defineAlphabet(char* string) {

	int* alphabet = calloc(256, sizeof * alphabet);

	if (alphabet) {

		for (unsigned int i = 0; i < strlen(string); i++) {
			int num = (unsigned)string[i];
			alphabet[num] = 1;
		}

		return alphabet;
	}
	else {
		fprintf(stderr, "Allocation of alphabet in defineAlphabet() failed.\n");
		freeAll();
		exit(0);
	}
}

/*
* The method adds the number of occurences of upper_case and lower_case of the same character
* and stores it in one variable for each character.
*
* int* components - an array that holds the number of different occurences of each character of the defined alphabet
*
* return - an array containing the total number of each character in the input, ignoring case
*/

int* addElements(const int* components) {

	//for molecular formulas only A-Z is used, 26 characters in total
	int* elementCount = calloc(26, sizeof(elementCount));

	if (elementCount) {

		for (unsigned int i = 0; i < strlen(ELEMENTS); i++) {
			// (int) 'A' == 65
			elementCount[ELEMENTS[i] - 65] = components[ORGANIC_SUBSET[2 * i]] + components[ORGANIC_SUBSET[(2 * i) + 1]];
		}

		*pAtoms = 0;
		//count the number of non-H atoms
		for (unsigned int i = 0; i < 26; i++) {

			if (i + 65 != 'H')
				numberOfAtoms += elementCount[i];
		}

		return elementCount;
	}
	else {
		fprintf(stderr, "Allocation of elementCount in addElements() failed.\n");
		freeAll();
		exit(0);
	}
}