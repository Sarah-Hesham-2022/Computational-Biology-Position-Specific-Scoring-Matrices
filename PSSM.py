'''
PSSM (Position Specific Scoring Matrix)

1. The multiple aligned DNA sequences is produced in 2 ways:

a. Reading the prepared multiple aligned DNA sequences from file called PSSMData.txt 
where first line is for the values of t sequences, and n nucleotides.

b. Generate t sequences where one sequence of length n is formed of random order of nucleotides.

2. Print the multiple aligned DNA sequences.

3. Apply PSSM on the multiple aligned DNA sequences.

4. Print the PSSM matrix

5. Let user enter a new sequence of length n, then print the probability of new entered sequence 
joining the rest of the multiple aligned DNA sequences.
'''

import random
import math

def randomSeqsGenterator(t,n):

     seqs = []
     for _ in range(t):
        seq = ""
        for _ in range(n): 
           seq += random.choice('ACGT')
        seqs.append(seq)

     return seqs

def prepareSequences():

    print("\nPlease Enter your choice for choosing multiple aligned DNA sequences:\n")
    print("1-Reading the prepared multiple aligned DNA sequences from file called PSSMData.txt where first line is for the values of t sequences, and n nucleotides.")
    print("2-Generate t sequences where one sequence of length n is formed of random order of nucleotides.\n")
    choice = 0
    t = 0
    n = 0
    seqs = []

    try:
       choice = int(input("Enter: \n"))
    except ValueError:
        print("Error, you are only allowed to enter two numbers 1 or 2.")
        return False

    if(choice < 1 or choice > 2):
        print("Error, you are only allowed to enter two numbers 1 or 2.")
        return False

    elif(choice == 1):
        try:
           myfile = open("PSSMData.txt")
           line = myfile.readline()
           line = line.split(" ")
           t = line[0]
           n = line[1]
           seqs = []
           line = myfile.readline()
           while(line):
              seqs.append(line.rstrip())
              line = myfile.readline() 
              
           t = len(seqs)
           n = len(seqs[0])

           for i in seqs:
              temp = len(i)
              if(temp != n):
                 print("Error: all sequences must be of same length, no insertions or deletions allowed.")
                 return False   
             
           return seqs

        except IndexError:
           print("Error, wrong file format.")
           return False
        except FileNotFoundError:
           print("Error, wrong file format.")
           return False
    else:
        try:
           t = int(input("Please, enter the t number for your sequences : \n"))
           n = int(input("Please, enter the n length number of all sequences : \n"))
           seqs = randomSeqsGenterator(t, n)
        except ValueError:
            print("Error, you are only allowed to enter numbers integers.")
            return False

        return seqs

def PSSM(seqs):

    t = len(seqs)
    n = len(seqs[0])
    array = []
    mergedSequences = ''.join(seqs)
    totalA = (mergedSequences.count("A") + mergedSequences.count("a")) / len(mergedSequences)
    totalC = (mergedSequences.count("C") + mergedSequences.count("c")) / len(mergedSequences)
    totalG = (mergedSequences.count("G") + mergedSequences.count("g")) / len(mergedSequences)
    totalT = (mergedSequences.count("T") + mergedSequences.count("t")) / len(mergedSequences)

    for i in range(4):
        temp = []
        for j in range(n):
            temp.append(0)
        array.append(temp)

    for i in range(n):
        A = 0.0
        C = 0.0 
        G = 0.0
        T = 0.0
        for j in range(t):
           if(seqs[j][i] == "a" or seqs[j][i] == "A"):
               A += 1
           elif(seqs[j][i] == "c" or seqs[j][i] == "C"):
               C += 1         
           elif(seqs[j][i] == "g" or seqs[j][i] == "G"):
               G += 1
           elif(seqs[j][i] == "t" or seqs[j][i] == "T"):
               T += 1

        array[0][i] = A/t/totalA
        array[1][i] = C/t/totalC
        array[2][i] = G/t/totalG
        array[3][i] = T/t/totalT

        if(array[0][i] != 0):
            array[0][i] = round(math.log(array[0][i], 2),2)
        if(array[1][i] != 0):
            array[1][i] = round(math.log(array[1][i], 2),2)
        if(array[2][i] != 0):
            array[2][i] = round(math.log(array[2][i], 2),2)
        if(array[3][i] != 0):
            array[3][i] = round(math.log(array[3][i], 2),2)

    return array

def calcProb(matrix,seq):

        prob = 0.0

        if(len(seq) != len(matrix[0])):
            print("\nError, your enetered sequence must be of the same length of the PSSM sequences length = " + str(len(matrix[0])) + "\n")
            return "Error"

        for i in range(len(seq)):
            if(seq[i] == "A" or seq[i] == "a"):
                prob += matrix[0][i]
            if(seq[i] == "C" or seq[i] == "c"):
                prob += matrix[1][i]        
            if(seq[i] == "G" or seq[i] == "g"):
                prob += matrix[2][i]
            if(seq[i] == "T" or seq[i] == "t"):
                prob += matrix[3][i]

        return round(prob,2)

def main():

    seqs = prepareSequences()
    matrix  = []

    if(seqs):
        print("\nThe multiple aligned DNA sequences:\n")
        for i in seqs:
            print(i)
        matrix = PSSM(seqs)
        print("\nThe PSSM matrix:\n")

        print('A: ',matrix[0])
        print('C: ',matrix[1])
        print('G: ',matrix[2])
        print('T: ',matrix[3])

        seq = input("\nEnter the sequence you want to calculate its probability : \n")

        probability = calcProb(matrix,seq)
        if(probability != "Error"):
           print("\nThe probability of new entered sequence joining the rest of the multiple aligned DNA sequences is : \n")
           print(probability)
           if(probability >= 0):
              print("\nThis sequence is a similar/identical residue match.\n")
           else:
              print("\nThis sequence is a non-conserved residue match.\n")

#################################################################################################################################
# Main Goes Here #
main()