import numpy as np
import random

rnd = random.Random()

dct  = { 'A':1 , 'T':2 , 'G':3 , 'C':4 }
dct2 = { 1:'A' , 2:'T' , 3:'G' , 4:'C' }

file = open(input('Enter File Containing DNA Sequences to be aligned using FTLPSO\n'))

NS,LN = [int(x) for x in input('Enter Number of DNA Sequences & length to be taken into consideration\n').split()]

mat = [[j for j in file.readline().strip()[:LN]]for i in range(NS)]

def match(st1,st2,KT,err):
    val = 0
    if len(st1)<KT or len(st2)<KT :
        return False
    for i in range(len(st1)):
        if st1[i]!=st2[i]:
            val+=1
        if val > err:
            break
    else:
        return True
    return False

KT,err = [int(x) for x in input('Enter Value of (K,Δ) for k-tuple table generation').split()]

score = np.inf

#Code for K,N-Tuple Algorithm, finding cluster of more uniform sizes
for n in range(NS):

    table = []
    indxs = [0]*NS

    for i in range(LN-KT+1):
        res = []
        try:
            for j in range(NS):
                if j!=n:
                    for k in range(indxs[j],LN):
                        if match(mat[n][i:i+KT] , mat[j][k:k+KT] , KT , err):
                            break
                    else:
                        raise Exception
                    res.append(k)
        except:
            continue
        if len(res) == NS-1:
            res.insert(n,i)
            table.append(res+[KT])
            for i in range(NS):
                indxs[i] = res[i]+1

    #code for adding the enteries of BEGIN & END
    table.insert(0 , [0]*NS + [0] )
    table.append( [LN-1]*NS + [0] )

    #Code For Merging the possible enteries
    while True:
        table2 = table[:]
        i=0
        while i<len(table)-1:
            if all( (table[i][j]+table[i][NS]) >= table[i+1][j] for j in range(NS) ):
                incrmnt = max( table[i+1][j] - table[i][j] for j in range(NS) )
                #incrmnt = max( table[i+1][j] + table[i+1][NS] - table[i][j] - table[i][NS] for j in range(NS) )
                table[i][NS] += incrmnt
                del table[i+1]
            else:
                i+=1
        if table == table2:
            break

    if all( (i+table[-1][-1]) > LN for i in table[-1][:-1]) :
        dcr = min( i + table[-1][-1] - LN for i in table[-1][:-1] )
        table[-1][-1] -= dcr

    #Code for Finding Score of current table, & making it main_table, if it has lesser score than, previous score
    new_score = 0
    for i in range(len(table)-1):
        for j in range(len(table[i])):
            new_score += (table[i+1][j]-table[i][j]+table[i][NS])**2

    new_score = new_score**0.5

    print('\nSequence Taken = {}, Score = {:0.3f} & Table having that score is followed'.format(n+1,new_score))
    for i in table:
        print(i)

    if new_score < score:
        score = new_score
        indx = n
        main_table = table

print('\n\nTable Generated at Iteration No. {} with Score of {:0.3f}'.format(indx+1,score))
for i in table:
    print(i)


#UPTO here Fragmentation is done completely, code ahead is for TLPSO
class Particle:
    def __init__(self,G,P):
        self.G = G
        self.P = P

    def distance(self,best):
        total_gaps = len(self.G) + len(best.G)
        matched_gaps = 0
        #code for finding Matched Gaps
        pass
        return ( matched_gaps / total_gaps )


'''
NP,R = [int(x) for x in input('Enter values of Particles in each Swarm & No. of Swarms in 1st Layer').split()]
c1 = 0.5    # Cognitive acceleration coefficient
c2 = 0.5    # Social    acceleration coefficient
'''
