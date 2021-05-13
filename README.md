# CheMGen

The program is a molecule generator which takes a molecular formula and a file as input an creates all possible molecules and writes it to the specified file.
The central method is recursive, isomorphic graphs are not rejected.

To use the program, download all files, go to the repository on cmd an type "make generator" (a C compiler is required).
After this the program is available, it can be called by the command "generator <formula> <outpuFile>", not including <> for the correct call.
Different data is written to the file: the input formula, the number of non-H atoms, the order of the atom identifiers corresponding to the row and column, each graph as upper triangle matrix on one line.
