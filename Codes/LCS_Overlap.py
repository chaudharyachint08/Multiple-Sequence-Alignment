# Dynamic Programming implementation of LCS problem
import time
import IP
import DNA

def lcs(X , Y):
    # find the length of the strings
    m = len(X)
    n = len(Y)
 
    # declaring the array for storing the dp values
    L = [[0 for x in range(n+1)] for x in range(m+1)]
 
    """Following steps build L[m+1][n+1] in bottom up fashion
    Note: L[i][j] contains length of LCS of X[0..i-1]
    and Y[0..j-1]"""
    for i in range(m+1):
        for j in range(n+1):
            if i == 0 or j == 0 :
                L[i][j] = 0
            elif X[i-1] == Y[j-1]:
                L[i][j] = L[i-1][j-1]+1
            else:
                L[i][j] = max(L[i-1][j] , L[i][j-1])
 
    # L[m][n] contains the length of LCS of X[0..n-1] & Y[0..m-1]
    return L[m][n]
#end of function lcs1


def lcs2(X, Y, m, n):
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
    lcs[index] = "\0"
 
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
 
    #print ( 'LCS Found is',"".join(lcs) )
#end of function lcs2 
 
# Driver program to test the above function


print('Length\t|LCS|\t|CLCS| \t Fractional Overhead\t  LCS%\t\t  CLCS%')
for L in range(3,14):
    X = IP.X[:2**L]
    Y = IP.Y[:2**L]
    #X = DNA.st
    #Y = X[::-1]

    #X='0123456789'
    #Y='5678901234'

    #X=''.join(chr(i) for i in range(ord('a'),ord('z')+1))+''.join(chr(i) for i in range(ord('A'),ord('Z')+1))
    #Y=''.join(chr(i) for i in range(ord('A'),ord('Z')+1))+''.join(chr(i) for i in range(ord('a'),ord('z')+1))

    #print('X=',X,'\nY=',Y)

    init = time.clock()
    #print( "Length of LCS is ", lcs(X, Y) )
    r1 = lcs(X,Y)
    #print( 'Seconds Taken',time.clock() - init )
    t1 = time.clock() - init

    if len(X)<len(Y):
        X=X*2
    else:
        Y=Y*2

    #print(len(X),len(Y))
    init = time.clock()
    #print( "Length of Circular LCS is ", lcs(X, Y) )
    r2 = lcs(X,Y)
    #print( 'Seconds Taken',time.clock() - init )
    print('{0:4}\t{1:4}   \t{2:4}\t\t{3:.4f}\t\t{4:.4f}%\t{5:.4f}%'.format(2**L,r1,r2,(time.clock() - init)/t1,(100*r1/(2**L)),(100*r2/(2**L))))
