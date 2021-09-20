'''
This is a simple python script to count the number
of different canconical SMILES representations which
are present in the in the input file.
'''

import sys

fileName = input("Please insert file name: ")

try:
    file = open(fileName, "r")

except FileNotFoundError:
    print("No such file or directory found.")
    sys.exit(1)    

if fileName.endswith(".can"):

    dictio = {}

    for line in file:
        if line in dictio:
            dictio[line] +=1
        else:
            dictio[line] = 1
    
    file.close()

    print(f"The number of different graphs is {len(dictio)}.")
    input()

else:
    print("Only files with extension '.can' are accepted.")
    input()
