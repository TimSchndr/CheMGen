#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<ctype.h>
#include<time.h>
#include "generator.h"

/*preprocessor variables to evaluate the algorithm:
usePartition: 0 - no composition and partition is used
			  1 - composition and partition is used

useMorgan: 0 - do not use the Morgan algorithm
		   1 - use the Morgan algorithm

weights: 0 - use the Morgan algorithm without taking edge weights into account
		 1 - use Morgan algorithm with edge weights
*/

#define usePartition 1
#define useMorgan 1
#define weights 1

// X stands for Chlorine
// B stands for Bromine
#define ORGANIC_SUBSET "CcNnOoPpSsBbFfHhIiXx"
#define ELEMENTS "CNOPSBFHIX"

//general information
int atomValence[26] = { 0 };

//counts for different formula properties
int numberOfnonHAtoms;
int numberOfAtoms;
int numberOfBonds;

//further properties
char* atomName;
int* atomIndex;
int* currentComposition;

//number of correct graphs
int counter = 0;

//the graph adjacency matrix
int* adjacencyMatrix;

//the Morgan matrices
long** morganNumbers;
int morganLength = 20;
int morganIndex = 0;

FILE* outputFile;

/*
* Driver Code
*/
void main(int argc, char* argv[]) {
	clock_t start, end;
	start = clock();

	//if no file is specified, write to stdout
	if (argc == 2) {
		outputFile = stdout;
	}
	else if (argc == 3) {
		outputFile = fopen(argv[2], "w");
	}
	else {
		fprintf(stderr, "Usage: generator <forumla> <outputFile>\n");
		exit(0);
	}

	if (outputFile) {

		//get the components of the formula
		int* components = getComponents(argv[1]);

		//count number of different nonH-atoms
		int differentAtoms = countAtomTypes(components);

		//set valence of atoms
		setAtomValence();

		//get the number of DBE
		int numberOfDBE = calcDoubleBonds(components);

		//get the number of bonds
		numberOfBonds = numberOfnonHAtoms - 1 + numberOfDBE;

		//define list of atoms
		defineAtomName(components, differentAtoms);

		free(components);

		//allocation of adjacencyMatrix
		adjacencyMatrix = calloc(numberOfnonHAtoms * numberOfnonHAtoms, sizeof(adjacencyMatrix));

		if (!adjacencyMatrix) {
			fprintf(stderr, "Allocation for adjacencyMatrix failed.\n");
			exit(0);
		}

#if useMorgan

		morganNumbers = calloc(morganLength, sizeof(morganNumbers));

		if (!morganNumbers) {
			fprintf(stderr, "Allocation failed for morganNumbers to store Morgan numbers.\n");
			exit(0);
		}
#endif

		int numberOfH = numberOfAtoms - numberOfnonHAtoms;

#if usePartition

		int** partitionList = malloc(differentAtoms * sizeof(partitionList));

		if (partitionList == NULL) {
			//allocation failed
			fprintf(stderr, "Allocation failed for partitionList in main.\n");
			exit(0);
		}

		//allocate partition array for each type of atoms
		for (int i = 0; i < differentAtoms; i++) {
			partitionList[i] = allocatePartition(atomIndex[i + 1] - atomIndex[i]);
		}

		while (getNextComposition(numberOfH, differentAtoms)) {
			//look at all possible compositions

			if (checkComposition(differentAtoms)) {
				//the composition is correct with respect to valence

				//number of generators
				int generatorCounter = differentAtoms - 1;
				//the generator to perform on
				int currentGenerator = differentAtoms - 1;

				for (int i = 0; i <= generatorCounter; i++) {
					//set all partitions to zero
					setPartitionToZero(atomIndex[i + 1] - atomIndex[i], partitionList[i]);
				}

				for (int i = 0; i < generatorCounter; i++) {
					//for all generators except the last get first partition
					nextAtomPartition(currentComposition[i], atomIndex[i + 1] - atomIndex[i], partitionList[i], atomName[atomIndex[i]]);
				}

				while (currentGenerator + 1) {

					int currentGenerator = differentAtoms - 1;

					while (nextAtomPartition(currentComposition[currentGenerator], atomIndex[currentGenerator + 1] - atomIndex[currentGenerator], partitionList[currentGenerator], atomName[atomIndex[currentGenerator]])) {
						//for the last generator try all partitions

						makeBond(0, 1, 0, partitionList, differentAtoms);

						if (currentComposition[currentGenerator] < 2) {
							//for less than two Hydrogens there is no other partition
							break;
						}
					}

					if (currentGenerator >= 0) {
						//there is more than one generator
						currentGenerator--;

						while ((currentGenerator >= 0 && !nextAtomPartition(currentComposition[currentGenerator], atomIndex[currentGenerator + 1] - atomIndex[currentGenerator], partitionList[currentGenerator], atomName[atomIndex[currentGenerator]])) || currentComposition[currentGenerator] == 0) {
							//search for next generator with correct partition
							currentGenerator--;
						}

						if (currentGenerator < 0) {
							//all generators do not have another partition
							break;
						}

						//another partition has been generated, go to the next one
						currentGenerator++;

						for (currentGenerator; currentGenerator < generatorCounter; currentGenerator++) {
							//for all generators behind the one with a new partition except the last
							setPartitionToZero(atomIndex[currentGenerator + 1] - atomIndex[currentGenerator], partitionList[currentGenerator]);
							nextAtomPartition(currentComposition[currentGenerator], atomIndex[currentGenerator + 1] - atomIndex[currentGenerator], partitionList[currentGenerator], atomName[atomIndex[currentGenerator]]);
						}

						//set last generator to zero
						setPartitionToZero(atomIndex[generatorCounter] - atomIndex[generatorCounter - 1], partitionList[generatorCounter]);
					}
				}
			}
		}
#else
		makeBond_1(0, 1, 0, differentAtoms);

#endif
		end = clock();
		double cpu_time_needed = (double)(end - start) / CLOCKS_PER_SEC;

		if (!argv[2]) {
			printf("Reported %d graph/s to stdout in %f seconds.\n", counter, cpu_time_needed);
		}
		else {
			printf("Reported %d graph/s to %s in %f seconds.\n", counter, argv[2], cpu_time_needed);
		}

		if (outputFile != stdout) {
			//only if there was an output file close it
			fclose(outputFile);
		}
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
* B - Bromine
* X - Chlorine
*/
void setAtomValence() {

	// (int) 'A' == 65
	atomValence['B' - 65] = 1;
	atomValence['C' - 65] = 4;
	atomValence['F' - 65] = 1;
	atomValence['H' - 65] = 1;
	atomValence['I' - 65] = 1;
	atomValence['N' - 65] = 3;
	atomValence['O' - 65] = 2;
	atomValence['P' - 65] = 2;
	atomValence['S' - 65] = 2;
	atomValence['X' - 65] = 1;

}

/*
* A method that returns the valence of the given atom
*
* char atom - the atom abbreviation for which the valence is required
* return - the atoms' valence
*/
int getAtomValence(char atom) {

	return atomValence[atom - 65];

}

/*
* A method that allocates an array with the atom identifier for each atom in the formula
* and an array which holds the number of each atom type which is not Hydrogen.
*/
void defineAtomName(int* components, int differentAtoms) {

	// +1 for terminal letter
	atomName = calloc(numberOfnonHAtoms + 1, sizeof(atomName));
	//+1 for terminal index
	atomIndex = calloc(differentAtoms + 1, sizeof(atomIndex));

	if (atomName && atomIndex) {

		//the index to fill atomName at each position
		int nameIndex = 0;
		//the index to fill atomIndex with the index of a new atom type
		int nextAtomType = 0;

		//set first index, always zero
		atomIndex[nextAtomType] = nameIndex;

		for (unsigned int i = 0; i < strlen(ELEMENTS); i++) {
			//amount of one atom type is processed, set next one if available

			if (ELEMENTS[i] != 'H' && components[ELEMENTS[i] - 65] > 0) {
				//for all components not 'H' fill the array with the atom names

				//there is another atom type
				nextAtomType++;

				for (int j = 0; j < components[ELEMENTS[i] - 65]; j++) {

					atomName[nameIndex] = ELEMENTS[i];
					nameIndex++;
				}
				atomIndex[nextAtomType] = nameIndex;
			}
		}
	}
	else {
		fprintf(stderr, "Allocation failed in defineAtomName().");
		freeAll();
		exit(0);
	}
}

/*
* A Method that checks if the valence of all atoms in the matrix is correct with respect
* to the calculated composition and partition.
*
* return - an int telling if all atoms are correct
*/
int checkAtoms(int** partitionList, int differentAtoms) {

	//index to look at in adjacencyMatrix
	int index = 0;

	for (int i = 1; i <= differentAtoms; i++) {
		for (int j = 0; j < atomIndex[i] - atomIndex[i - 1]; j++) {

			if (getAtomValence(atomName[index]) - matrix(index, index) != partitionList[i - 1][j]) {
				//for currently observed atom valence is not correct
				return FALSE;
			}
			index++;
		}
	}
	return TRUE;
}

int checkAtoms_1(int differentAtoms) {

	int index = 0;
	for (int i = 0; i <= differentAtoms; i++) {

		if (getAtomValence(atomName[index]) < matrix(index, index)) {
			//for currently observed atom valence is not correct
			return FALSE;
		}
		index++;
	}
	return TRUE;
}

/*
* A method that fills the adjacencyMatrix with bonds recursively.
*
* int i - the row index
* int j - the column index
* int bondSum - the number of bonds already spent
* int** partition - the current hydrogen partition
* int differentAtoms - the number of different atoms
*/
void makeBond(int i, int j, int bondSum, int** partition, int differentAtoms) {

	if (numberOfBonds == bondSum) {
		// all bonds have been placed

		if (checkAtoms(partition, differentAtoms)) {
			// check if all atoms are ok

			if (checkConnectivity()) {
				//check if the graph is connected
#if useMorgan
				long* molecule = calculateMorganNumbers();

				if (!morganNumbersExist(molecule)) {

					saveMorganNumbers(molecule);
#endif
					//save the graph
					writeSDFformat();
					counter++;
#if useMorgan
				}
#endif
			}
		}
		return;
	}

	if ((i * numberOfnonHAtoms + j) >= (numberOfnonHAtoms * numberOfnonHAtoms)) {
		//end of the matrix is reached
		return;
	}

	//add a new bond between i and j
	for (int b = 0; b <= MIN_VALUE(getAtomValence(atomName[i]), getAtomValence(atomName[j])); b++) {

		if (bondSum + b > numberOfBonds) {
			break;
		}

		matrix(i, j) = b;
		matrix(j, i) = b;
		matrix(i, i) += b;
		matrix(j, j) += b;

		//gesamt - belegte = freie >= benötigt

		if (getAtomValence(atomName[i]) - matrix(i, i) >= getPartitionCount(i, partition) && getAtomValence(atomName[j]) - matrix(j, j) >= getPartitionCount(j, partition)) {
			//check if valence for both atoms is ok with respect to partition

			if (j < numberOfnonHAtoms - 1) {
				//same row, next column
				makeBond(i, j + 1, bondSum + b, partition, differentAtoms);
			}
			else {
				//start with the next row
				makeBond(i + 1, i + 2, bondSum + b, partition, differentAtoms);
			}
		}
		else {
			matrix(i, i) -= b;
			matrix(j, j) -= b;
			break;

		}
		matrix(i, i) -= b;
		matrix(j, j) -= b;
	}

	matrix(i, j) = 0;
	matrix(j, i) = 0;
}

/*
* A method that fills the adjacency matrix but withou respect to a partition of the number of hydrogens.
*
* mainly run time checking option
*/
void makeBond_1(int i, int j, int bondSum, int differentAtoms) {

	if (numberOfBonds == bondSum) {
		// all bonds have been placed

		if (checkAtoms_1(differentAtoms)) {
			// check if all atoms are ok

			if (checkConnectivity()) {
				//check if the graph is connected
#if useMorgan
				long* molecule = calculateMorganNumbers();

				if (!morganNumbersExist(molecule)) {

					saveMorganNumbers(molecule);
#endif
					//save the graph
					writeSDFformat();
					counter++;
#if useMorgan
				}
#endif
			}
		}
		return;
	}

	if ((i * numberOfnonHAtoms + j) >= (numberOfnonHAtoms * numberOfnonHAtoms)) {
		//end of the matrix is reached
		return;
	}

	//add a new bond between i and j
	for (int b = 0; b <= MIN_VALUE(getAtomValence(atomName[i]), getAtomValence(atomName[j])); b++) {

		if (bondSum + b > numberOfBonds) {
			break;
		}

		matrix(i, j) = b;
		matrix(j, i) = b;
		matrix(i, i) += b;
		matrix(j, j) += b;

		if (matrix(i, i) <= getAtomValence(atomName[i]) && matrix(j, j) <= getAtomValence(atomName[j])) {
			//check if valence for both atoms is ok

			if (j < numberOfnonHAtoms - 1) {
				//same row, next column
				makeBond_1(i, j + 1, bondSum + b, differentAtoms);
			}
			else {
				//start with the next row
				makeBond_1(i + 1, i + 2, bondSum + b, differentAtoms);
			}
		}
		else {
			matrix(i, i) -= b;
			matrix(j, j) -= b;
			break;

		}
		matrix(i, i) -= b;
		matrix(j, j) -= b;
	}

	matrix(i, j) = 0;
	matrix(j, i) = 0;
}

/*
* A method that does a DFS on an adjacency matrix for given start node.
*
* int* visitedNodes - an array holding a boolean to tell if a node has been visited before
*/
void dfs(int start, int* visitedNodes) {

	visitedNodes[start] = 1;

	//only if the start node has degree > 0 the Graph can be connected
	if (matrix(start, start)) {

		//search for the next adjacent node to start from (not on main diagonal)
		for (int i = 0; i < numberOfnonHAtoms; i++) {
			//if vertex i has not been visited and there is an edge from start to i, do dfs
			if (!visitedNodes[i] && matrix(start, i)) {
				dfs(i, visitedNodes);
			}
		}
	}
}

/*
* A method that checks if the underlying graph is connected.
*/
int checkConnectivity() {

	int* visitedNodes = calloc(numberOfnonHAtoms, sizeof(visitedNodes));

	if (visitedNodes) {

		//start the DFS with node 0
		dfs(0, visitedNodes);

		//check if all nodes have been visited
		for (int i = 0; i < numberOfnonHAtoms; i++) {
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
A method that finds the smallest cycle in an undirected graph, ignoring edge weights.

return - the length of the shortest cycle in the graph
*/
int shortestCycle() {

	int maxDist = 100;

	//set maximum cycle length
	int length = maxDist;

	//array to store distances
	int* distance = calloc(numberOfnonHAtoms, sizeof(distance));

	if (!distance) {
		fprintf(stderr, "Allocation of distance array failed in method shortestCycle.");
		freeAll();
		exit(0);
	}

	//notice parent for each node
	int* parent = calloc(numberOfnonHAtoms, sizeof(parent));

	if (!parent) {
		fprintf(stderr, "Allocation of parent array failed in method shortestCycle.");
		freeAll();
		exit(0);
	}

	//representation of a queue
	int* queue = calloc(2 * numberOfnonHAtoms, sizeof(queue));

	if (!queue) {
		fprintf(stderr, "Allocation of queue array failed in method shortestCycle.");
		freeAll();
		exit(0);
	}

	//do the procedure for every vertex
	for (int v = 0; v < numberOfnonHAtoms; v++) {

		//read-write indices for queue
		int insertIndex = 0;
		int readIndex = 0;

		for (int i = 0; i < numberOfnonHAtoms; i++) {
			//set distances maximum
			distance[i] = maxDist;
			//set parent to imaginary
			parent[i] = -1;
		}

		//push current node to queue
		queue[insertIndex] = v;
		insertIndex++;

		//set distance of v to 0
		distance[v] = 0;

		//queue is not empty
		while (insertIndex != readIndex) {

			//read entry from queue
			int x = queue[readIndex];
			readIndex++;

			for (int j = 0; j < numberOfnonHAtoms; j++) {

				if ((j != x) && matrix(x, j)) {
					//get all neighbours of current vertex x

					if (distance[j] == maxDist) {
						//j was not found yet

						//set distance to distance from x + 1
						distance[j] = distance[x] + 1;

						parent[j] = x;

						queue[insertIndex] = j;
						insertIndex++;
					}
					else if (parent[x] != j && parent[j] != x) {
						//closed a cycle, x and j are not neighbours

						//report length of the cycle
						length = MIN_VALUE(length, (distance[x] + distance[j] + 1));
					}
				}
			}
		}
	}

	free(queue);
	free(parent);
	free(distance);

	if (length == maxDist) {
		//no cycle in graph
		return -1;
	}
	else {
		//length of shortest cycle
		return length;
	}

}

/*
* A method that creates the next weak composiiton on basis of the last one.
*
* int n - the integer to decompose
* int k - the number of parts for the composition
*
* return - int telling True or False
*/
int getNextComposition(int n, int k) {

	if (k <= 0) {
		//there is no box
		fprintf(stdout, "Without a box, no composition.\n");
		return FALSE;
	}
	else {
		if (currentComposition == NULL) {
			//allocate array
			currentComposition = calloc(sizeof(currentComposition), k);

			if (currentComposition == NULL) {
				//check allocation
				fprintf(stdout, "Allocation failed in yieldComp().\n");
				exit(0);
			}

			//set initial composition
			currentComposition[k - 1] = n;
			return TRUE;
		}
		else {
			//array is allocated, get next Composition

			if (k == 1) {
				//only one composition, has been generated
				return FALSE;
			}

			if (currentComposition[0] == n) {
				//last composiiton was the final one, end
				return FALSE;
			}

			if (currentComposition[k - 1] != 0) {
				//remainder not zero, generate next composition
				currentComposition[k - 1]--;
				currentComposition[k - 2]++;
				return TRUE;
			}
			else {
				//start search for highest position in front of remainder
				for (int i = k - 2; i >= 0; i--) {
					if (currentComposition[i] != 0) {
						currentComposition[i - 1]++;
						currentComposition[i] = 0;
						break;
					}
				}
				//count spent values
				int sum = 0;
				for (int i = 0; i < k - 1; i++) {
					sum += currentComposition[i];
				}

				//reset remainder
				currentComposition[k - 1] = n - sum;
				return TRUE;
			}
		}
	}
}

/*
* A method that checks if the the number of atoms for each type
* can have the amount of valent hydrogens.
*
* return - an int telling True or False
*/
int checkComposition(int differentAtoms) {

	if (differentAtoms == 1) {
		if (atomName[0] == 'C') {
			if (currentComposition[0] == 4) {
				return TRUE;
			}
		}
		else if (atomName[0] == 'N') {
			if (currentComposition[0] == 3) {
				return TRUE;
			}
		}
		else {
			return FALSE;
		}
	}

	for (int i = 0; i < differentAtoms; i++) {

		if (currentComposition[i] > ((atomIndex[i + 1] - atomIndex[i]) * (getAtomValence(atomName[atomIndex[i]]) - 1) - (atomIndex[i + 1] - atomIndex[i] - 2))) {
			return FALSE;
		}
	}
	return TRUE;
}

/*
* A method that allocates an array to hold an integer partition.
*
* int partitionLength - the number of atoms of the type
*
* return - the allocated array
*/
int* allocatePartition(int partitionLength) {

	int* partition = calloc(partitionLength, sizeof(partition));

	if (partition == NULL) {
		//check allocation
		fprintf(stdout, "Allocation failed in allocatePartition().\n");
		exit(0);
	}

	return partition;
}

/*
* A method that replaces the current partition with zeros.
*
* int length - the length of the partition
* int* partititon - the array that holds the partition
*/
void setPartitionToZero(int length, int* partition) {

	if (length < 1) {
		fprintf(stderr, "Length of partition must be greater than 0.\n");
		exit(0);
	}

	for (int i = 0; i < length; i++) {
		partition[i] = 0;
	}

}

/*
*	A method that calculates a new partition on basis of the given on.
*
* int number - the integer to be partitioned
* int* array - the array that holds the partition
*
*  adopted from https://www.sanfoundry.com/c-program-perform-partition-integer-possible-ways/ (04.07.2021)
*/
int nextPartition(int number, int* array) {

	int currentIndex = number - 1;

	//find last index in currentPartition with non-0 value
	for (currentIndex; array[currentIndex] == 0; currentIndex--);

	int remValue = 0;

	//note all 1s and set them to 0
	while (currentIndex >= 0 && array[currentIndex] == 1) {
		remValue += array[currentIndex];
		array[currentIndex] = 0;
		currentIndex--;
	}

	if (currentIndex < 0) {
		//last partition, only consisting of 1s
		return FALSE;
	}

	//decrease number at found index
	array[currentIndex]--;
	remValue++;

	while (remValue > array[currentIndex]) {
		array[currentIndex + 1] = array[currentIndex];
		remValue = remValue - array[currentIndex];
		currentIndex++;
	}

	//set remValue to next position and increment index
	array[currentIndex + 1] = remValue;
	currentIndex++;

	return TRUE;
}

/*
* A method that calculates the next integer partition of a number
* with respect to the current partition.
*
* int numberOfH - the number of Hydrogens to be partitioned
* int length - the length of the array
* int* partition - the array that holds the partition
* char atomType - the atom type for the partition
*
* return - an int telling wheter the current partition is ok
*/
int nextAtomPartition(int numberOfH, int arrayLength, int* partition, char atomType) {

	if (arrayLength < 1) {
		fprintf(stderr, "int arrayLength must be greater than 0.\n");
		exit(0);
	}

	if (partition[0] == 0) {
		//no partition set

		if (numberOfnonHAtoms == 1) {
			//if only one atom exists
			if (atomName[0] == 'C') {
				partition[0] = 4;
				return TRUE;
			}
			else if (atomName[0] == 'N') {
				partition[0] = 3;
				return TRUE;
			}
		}

		int number = numberOfH;
		int valence = getAtomValence(atomType) - 1;

		//set the first partition
		for (int i = 0; i < arrayLength; i++) {
			if (number >= valence) {
				partition[i] = valence;
				number -= valence;
			}
			else {
				partition[i] = number;
				number = 0;
			}
		}

		return TRUE;
	}

	else if (partition[0] == 1) {
		//last partition
		return FALSE;
	}

	else {

		//array to hold the complete partition of numberOfH
		int* wholePartition = calloc(numberOfH, sizeof(partition));

		if (wholePartition == NULL) {
			fprintf(stderr, "Allocation failed in nextAtomPartition().\n");
			exit(0);
		}

		//copy the current partition into arrayLength parts into wholePartition
		memcpy(wholePartition, partition, arrayLength * sizeof(wholePartition));

		while (nextPartition(numberOfH, wholePartition)) {
			//calculate next partitions of numberOfH

			if (wholePartition[arrayLength] == 0 || arrayLength == numberOfH) {
				//there is a partition of numberOfH into arrayLength parts

				//copy wholePartition into the array partition with arrayLength parts
				memcpy(partition, wholePartition, arrayLength * sizeof(wholePartition));

				free(wholePartition);
				return TRUE;
			}
		}
		free(wholePartition);
		return FALSE;
	}
}

/*
* A method that checks if the current partition can be applied to the atomName array.
*
* return - an int that tells whether the partition can be applied
*/
int checkAtomPartition(int number, int* partition) {

	int lastIndex = number - 1;

	//search last index of partition
	for (lastIndex; partition[lastIndex] == 0; lastIndex--);

	if (lastIndex + 1 > numberOfnonHAtoms) {
		//more numbers in partition than nonH-atoms in molecule
		return FALSE;
	}

	for (int i = 0; i < numberOfnonHAtoms; i++) {
		if (getAtomValence(atomName[i]) > partition[i]) {
			//the atom must be connected and saturated
			continue;
		}
		else if (numberOfnonHAtoms == 1 && partition[i] == getAtomValence(atomName[i])) {
			//only one atom
			return TRUE;
		}
		else {
			return FALSE;
		}
	}

	//all atoms can have the number of H-atoms
	return TRUE;
}

/*
* A method that returns the number of hydrogens assigned to the atom of interest.
*
* int index - the atom index of interest
* int** the partition list
* int differentAtoms - the number of different atoms
*
* return - the number of hydrogens assigned to atom index
*/
int getPartitionCount(int index, int** partition) {

	int pos = 1;

	while (atomIndex[pos] <= index) {
		pos++;
	}
	pos--;

	return partition[pos][index - atomIndex[pos]];
}

/*
* This method writes the underlying chemical graph
* to the outputFile in Molfile format.
* To meet the format, the hydrogens are enumerated and
* attached to one specific atom.
*/
void writeSDFformat() {

	//calculate the number of internal connections
	int connections = 0;

	for (int i = 0; i < numberOfnonHAtoms; i++) {
		for (int j = i + 1; j < numberOfnonHAtoms; j++) {
			if (matrix(i, j)) {
				connections++;
			}
		}
	}

	fprintf(outputFile, "\nCheMGen 1.0\nMolecule %d\n", counter + 1);

	//Count Line Block
	fprintf(outputFile, "%3d%3d  0  0  0  0            999 V2000\n", numberOfAtoms, connections + numberOfAtoms - numberOfnonHAtoms);

	//Atom Block
	for (int i = 0; i < numberOfnonHAtoms; i++) {
		//write non-H atoms to file
		fprintf(outputFile, "    0.0000    0.0000    0.0000 %c  0  0  0  0  0  0  0  0  0  0  0  0\n", atomName[i]);
	}

	for (int i = numberOfnonHAtoms + 1; i <= numberOfAtoms; i++) {
		//write all H atoms to file
		fprintf(outputFile, "    0.0000    0.0000    0.0000 %s  0  0  0  0  0  0  0  0  0  0  0  0\n", "H");
	}

	//Bond Block for non-H atoms
	for (int i = 0; i < numberOfnonHAtoms; i++) {
		for (int j = i + 1; j < numberOfnonHAtoms; j++) {

			//if there is a bond, report to file
			if (matrix(i, j)) {
				fprintf(outputFile, "%3d%3d%3d  0  0  0  0\n", i + 1, j + 1, matrix(i, j));
			}
		}
	}

	//Bond Block for H atoms

	//represents the current H-atom number
	int hCounter = numberOfnonHAtoms + 1;

	for (int i = 0; i < numberOfnonHAtoms; i++) {
		for (int b = matrix(i, i); b < getAtomValence(atomName[i]); b++) {
			//parameters: non-H atom, H-atom, bond multiplicity
			fprintf(outputFile, "%3d%3d%3d  0  0  0  0\n", i + 1, hCounter, 1);

			hCounter++;
		}
	}
	//End of File
	fprintf(outputFile, "M  END\n$$$$\n");
}

/*
* Free all global variables.
* This method is called when an error occured during run.
*/
int freeAll() {

	if (atomName) {
		free(atomName);
	}
	if (adjacencyMatrix) {
		free(adjacencyMatrix);
	}

	return 1;
}

/*
* A method that calculates the morgan numbers for the current adjacency matrix and returns it.
*
* return - a pointer to long array containing the Morgan numbers
*/
long* calculateMorganNumbers() {

	//arrays to hold old and new morgan numbers
	long* currentMorganNumbers = calloc(numberOfnonHAtoms, sizeof(currentMorganNumbers));
	long* newMorganNumbers = calloc(numberOfnonHAtoms, sizeof(newMorganNumbers));

	if (!currentMorganNumbers || !newMorganNumbers) {
		fprintf(stderr, "Allocation failed in calculateMorganNumbers.\n");
		exit(0);
	}

	//copy array to have an array to sort
	long* copyArray = calloc(numberOfnonHAtoms, sizeof(copyArray));

	if (!copyArray) {
		fprintf(stderr, "Allocation failed for copyArray in calculateMorganNumbers().\n");
		exit(0);
	}

	//ints to count the number of different digits in the arrays
	int oldCount = 0;
	int newCount = 0;

	//set the first Morgan number for each atom (number of neighbours)
	for (int i = 0; i < numberOfnonHAtoms; i++) {
		for (int j = 0; j < numberOfnonHAtoms; j++) {
			if (i != j && matrix(i, j)) {
				currentMorganNumbers[i] += 1;
			}
		}
	}

	//count number of different values
	copyArray = memcpy(copyArray, currentMorganNumbers, numberOfnonHAtoms * sizeof(currentMorganNumbers));
	quickSort(copyArray, copyArray + numberOfnonHAtoms - 1);
	oldCount = countDigits(copyArray, numberOfnonHAtoms);

	for (int i = 0; i < numberOfnonHAtoms; i++) {
		newMorganNumbers[i] = 0;
		for (int j = 0; j < numberOfnonHAtoms; j++) {
			if (i != j && matrix(i, j) > 0) {
				//check if two atoms are adjacent

#if weights
				newMorganNumbers[i] += matrix(i, j) * currentMorganNumbers[j];
#else
				newMorganNumbers[i] += currentMorganNumbers[j];
#endif			
			}
		}
	}

	copyArray = memcpy(copyArray, newMorganNumbers, numberOfnonHAtoms * sizeof(newMorganNumbers));
	quickSort(copyArray, copyArray + numberOfnonHAtoms - 1);
	newCount = countDigits(copyArray, numberOfnonHAtoms);

	while (newCount > oldCount) {

		//set new morgan numbers to old
		memcpy(currentMorganNumbers, newMorganNumbers, numberOfnonHAtoms * sizeof(currentMorganNumbers));

		oldCount = newCount;

		for (int i = 0; i < numberOfnonHAtoms; i++) {
			newMorganNumbers[i] = 0;
			for (int j = 0; j < numberOfnonHAtoms; j++) {
				if (i != j && matrix(i, j) > 0) {
					//check if two atoms are adjacent
					newMorganNumbers[i] += currentMorganNumbers[j];
				}
			}
		}

		copyArray = memcpy(copyArray, newMorganNumbers, numberOfnonHAtoms * sizeof(newMorganNumbers));
		//quickSort(copyArray, copyArray + numberOfnonHAtoms - 1);
		newCount = countDigits(copyArray, numberOfnonHAtoms);
	}

	memcpy(currentMorganNumbers, newMorganNumbers, numberOfnonHAtoms * sizeof(currentMorganNumbers));

	free(newMorganNumbers);
	free(copyArray);
	//quickSort(currentMorganNumbers, currentMorganNumbers + numberOfnonHAtoms - 1);
	return currentMorganNumbers;
}

/*
* A method that checks if the Morgan numbers of the given molecule already exist.
*
* long* molecule - the currently calculated Morgan numbers
*
* return - an int if the Morgan numbers already exist
*/
int morganNumbersExist(long* molecule) {

	for (int i = 0; i < morganIndex; i++) {
		//look at all Morgan numbers
		int j = 0;
		for (j = 0; j < numberOfnonHAtoms; j++) {
			//compare each position in the arrays
			if (molecule[j] != morganNumbers[i][j]) {
				break;
			}
		}

		if (j == numberOfnonHAtoms) {
			//second loop completed without break, Morgan numbers exist
			return TRUE;
		}
	}
	return FALSE;
}

/*
* A method that stores a new Morgan matrix into the pointer array morganNumbers at the right position.
*
* int index - the index to put the Morgan numbers in
* long* molecule - the Morgan numbers
*/
void saveMorganNumbers(long* molecule) {

	if (morganIndex == morganLength) {
		//the array is full, must be reallocated

		morganLength = morganLength * 2;
		long** temp = realloc(morganNumbers, morganLength * sizeof(temp));

		if (temp) {
			morganNumbers = temp;
		}
		else {
			fprintf(stderr, "Reallocation failed for morganNumbers in saveMorganNumbers().\n");
			exit(0);
		}
	}

	morganNumbers[morganIndex] = molecule;

#if debug
	for (int i = 0; i < morganIndex; i++) {
		for (int j = 0; j < numberOfnonHAtoms; j++) {
			printf("%d ", morganNumbers[i][j]);
		}
		printf("\n");
	}
#endif

	morganIndex++;
	return;
}

/*
* A method that takes an increasingly sorted array and counts the number of different digits.
*
* int* array - the sorted array in increasing order
* int length - the length of the array
*
* return - the number of different digits
*/
int countDigits(long* array, int length) {

	int count = 1;

	for (int i = 1; i < length; i++) {
		if (array[i] == array[i - 1]) {
			continue;
		}
		else {
			count++;
		}
	}

	return count;
}

//Methods for the formula reader

/*
* A method that reads in a String from console (molecular formula) and
* counts the number of occurences of each element.
*
* return - an array containing the count for each letter in the input
*/
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

					memcpy(num, &formula[i], len + 1 * sizeof(num));

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
int calcDoubleBonds(const int* elementCounts) {

	//total number of halogens in the input
	int halogenSum = 0;
	char halogens[5] = "BFIX";

	for (unsigned int i = 0; i < strlen(halogens); i++) {
		halogenSum += elementCounts[halogens[i] - 65];
	}

	//total number of two-valent atoms
	int twoVal = 0;
	char twoValAtoms[4] = "OPS";

	for (unsigned int i = 0; i < strlen(twoValAtoms); i++) {
		twoVal += elementCounts[twoValAtoms[i] - 65];
	}

	//check if only one structural atom (C or N) exists
	if (elementCounts['C' - 65] == 1 && elementCounts['N' - 65] == 0 && !twoVal) {

		//check if the atom is saturated
		if (halogenSum + elementCounts['H' - 65] != getAtomValence('C')) {

			fprintf(stderr, "The carbon atom was not saturated.\n");
			freeAll();
			exit(0);
		}
	}
	else if (elementCounts['C' - 65] == 0 && elementCounts['N' - 65] == 1 && !twoVal) {

		//check if the atom is saturated
		if (halogenSum + elementCounts['H' - 65] != getAtomValence('N')) {

			fprintf(stderr, "The nitrogen atom was not saturated.\n");
			freeAll();
			exit(0);
		}
	}

	if (numberOfAtoms == 2 && elementCounts['C' - 65] == 2) {
		//check if only two carbon atoms exist
		fprintf(stderr, "No valid structures possible for this input.\n");
		exit(0);
	}

	if (elementCounts['C' - 65] == 4 && elementCounts['H' - 65] == 2) {
		fprintf(stderr, "This formula is not processable at the moment.\n");
		exit(0);
	}

	//the number of H, N and halogens combined according to DBE formula
	double negativeTerm = 0;
	negativeTerm = 0.5 * ((double)elementCounts['H' - 65] - (double)elementCounts['N' - 65] + (double)halogenSum);
	int dbe = -1;

	if (negativeTerm == (int)negativeTerm) {
		//check if negativeTerm is an int
		dbe = elementCounts['C' - 65] + 1 - negativeTerm;
	}

	if (dbe >= 0) {
		return dbe;
	}
	else {
		fprintf(stderr, "The input was not correct with respect to DBE calculation.\n");
		freeAll();
		exit(0);
	}
}

/*
* A method that counts the number of different atom types except Hydrogen.
*
* int* components - an array that holds the counts for each atom
*
* return - the number of different atom types
*/
int countAtomTypes(int* components) {

	int atomTypes = 0;

	for (int i = 0; i < 26; i++) {

		if (components[i] && i + 65 != 'H') {
			atomTypes++;
		}
	}

	return atomTypes;
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

		numberOfAtoms = 0;
		numberOfnonHAtoms = 0;

		//count the number of non-H atoms
		for (unsigned int i = 0; i < 26; i++) {

			numberOfAtoms += elementCount[i];

			if (i + 65 != 'H') {
				numberOfnonHAtoms += elementCount[i];
			}

		}
		return elementCount;
	}
	else {
		fprintf(stderr, "Allocation of elementCount in addElements() failed.\n");
		freeAll();
		exit(0);
	}
}

//QuickSort method

/*
* An implementation of the Quicksort sorting algorithm.
* adopted from https://openbook.rheinwerk-verlag.de/c_von_a_bis_z/022_c_algorithmen_003.htm (11.08.21
*/
void quickSort(long* left, long* right) {
	long* ptr1 = left;
	long* ptr2 = right;
	long w, mid;

	mid = *(left + (right - left >> 1));

	do {
		while (*ptr1 < mid) {
			ptr1++;
		}
		while (*ptr2 > mid) {
			ptr2--;
		}

		if (ptr1 > ptr2) {
			break;
		}

		w = *ptr1;
		*ptr1 = *ptr2;
		*ptr2 = w;

	} while (++ptr1 <= --ptr2);
	if (left < ptr2) {
		quickSort(left, ptr2);
	}
	if (ptr1 < right) {
		quickSort(ptr1, right);
	}
}