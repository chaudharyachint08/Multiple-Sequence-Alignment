import numpy as np
import random

rnd = random.Random()

dct  = { 'A':1 , 'T':2 , 'G':3 , 'C':4 }
dct2 = { 1:'A' , 2:'T' , 3:'G' , 4:'C' }


def lcs2(X, Y):

    m = len(X)
    n = len(Y)
    
    L = [[0 for x in range(n+1)] for x in range(m+1)]
 
    # Following steps build L[m+1][n+1] in bottom up fashion. Note
    # that L[i][j] contains length of LCS of X[0..i-1] and Y[0..j-1] 
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0:
                L[i][j] = 0
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1] + 1
            else:
                L[i][j] = max(L[i-1][j], L[i][j-1])
 
    # Following code is used to print LCS
    index = L[m][n]
 
    # Create a character array to store the lcs string
    lcs = [""] * (index+1)
    lcs[index] = ""
 
    # Start from the right-most-bottom-most corner and
    # one by one store characters in lcs[]
    i = m
    j = n
    while i > 0 and j > 0:
 
        # If current character in X[] and Y are same, then
        # current character is part of LCS
        if X[i-1] == Y[j-1]:
            lcs[index-1] = X[i-1]
            i-=1
            j-=1
            index-=1
 
        # If not same, then find the larger of two and
        # go in the direction of larger value
        elif L[i-1][j] > L[i][j-1]:
            i-=1
        else:
            j-=1
 
    return "".join(lcs)
#end of function lcs2 


file = open(input('Enter File Containing DNA Sequences to be aligned using Classical LCS\n'))

NS,LN = [int(x) for x in input('Enter Number of DNA Sequences & length to be taken into consideration\n').split()]

mat = [file.readline().strip()[:LN] for i in range(NS)]
mat2 = mat[:]
mat3 = mat[:]


#LOOP FOR GREEDY FROM BEGINNING BASED CONSENSUS EVALUTION
count=0
while len(mat3)!=1:
    count+=len(mat3)
    mat3.append(lcs2(mat3.pop(0),mat3.pop(0)))
consensus = mat3[0]
print('LCS Based Greedy Approach with length',len(consensus),'\n',consensus)

'''
#LOOP FOR GREEDY FROM LAST BASED CONSENSUS EVALUTION
count=0
while len(mat2)!=1:
    count+=len(mat2)
    mat2.append(lcs2(mat2.pop(),mat2.pop()))

consensus = mat2[0]

print('\nLCS Based Greedy from LAST Approach with length',len(consensus),'\n',consensus)
'''

#LOOP FOR TOURNAMENT BASED CONSENSUS EVALUTION
count=0
while len(mat)!=1:
    count+=len(mat)
    mat = [lcs2(mat[i],mat[i+1]) for i in range( 0 , len(mat)-( 1 if len(mat)%2==1 else 0 ) , 2 ) ] + ( [mat[-1]] if len(mat)%2==1 else [] )

consensus = mat[0]

print('\nLCS Based Tournament Approach with length',len(consensus),'\n',consensus)
