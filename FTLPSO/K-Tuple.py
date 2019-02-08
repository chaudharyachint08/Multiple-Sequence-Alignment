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

#    print('\nSequence Taken = {}, Score = {:0.3f} & Table having that score is followed'.format(n+1,new_score))
#    for i in table:
#        print(i)

    if new_score < score:
        score = new_score
        indx = n
        main_table = table

print('\nTable Generated at Iteration No. {} with Score of {:0.3f}\n'.format(indx+1,score))
#for i in table:
#    print(i)

aligned = [[] for i in range(NS)]

for i in range(len(table)-1):

    seq1 = [ mat[j][ table[i][j] : table[i][j]+table[i][NS] ]  for j in range(NS) ]
    seq2 = [ mat[j][ table[i][j]+table[i][NS] : table[i+1][j] ]  for j in range(NS) ]
    
        
    ln = max( len(i) for i in seq2 )

    for i in seq2:
        i.extend(['-']*(ln-len(i)))

    for j in range(NS):
        temp = np.array(seq2[j])
        np.random.shuffle(temp)
        temp = [x for x in temp]
        aligned[j].append( ''.join(seq1[j]) )
        aligned[j].append( ''.join(temp) )


lst = [ mat[j][ table[-1][j] : table[-1][j]+table[-1][NS] ]  for j in range(NS) ]

for i in range(NS):
    aligned[i].append( ''.join(lst[i]) )


def show_aligned():
    for i in range(NS):
        print(aligned[i])


#show_aligned()

max_len = 0

for i in range(len(aligned)):
    aligned[i] = ''.join(aligned[i])
    val = len(aligned[i])
    if val > max_len:
        max_len = val

for i in range(len(aligned)):
    aligned[i] += '-'*(max_len-len(aligned[i]))

def print_msa():
    for i in aligned:
        print(i)


#print(aligned[0])
mul = (LN/len(aligned[0]))**1.5


def CS():

    score = 0

    for i in range(len(aligned[0])):
        su = { 'A':0 , 'T':0 , 'G':0 , 'C':0 , '-':0 }

        for j in aligned:
            su[j[i]] += 1

        m = max(su.values())
        score += (m*(1+m/LN))
    
    return score*mul


def SOP():
    '''Sum_Of_Pair Function'''

    score = 0

    nA,nB,nC = 0,0,0
    A,B,C    = 9,-5,-2

    for i in aligned:
        for j in aligned:
            if i is not j:
                for k in range(len(aligned[0])):
                    if i[k]==j[k]:
                        nA += 1
                    elif '-' in (i[k],j[k]):
                        nB += 1
                    else:
                        nC += 1

    score = A*nA - B*nB - C*nC

    return str(round(nA*mul))+'*A - '+str(round(nB*mul))+'*B - '+str(round(nC*mul))+'*C' , int(round(score*mul))

consensus = []
for i in range(len(aligned[0])):
    su = { 'A':0 , 'T':0 , 'G':0 , 'C':0 , '-':0 }

    for j in range(NS):
        su[ aligned[j][i] ] += 1

    su = list(su.items())
    su.sort(key=lambda x:x[1],reverse=True)

    if su[0][1] >= su[1][1]*(1/mul**(2/3)):
        if su[0][0] != '-':
            consensus.append(su[0][0])

consensus = ''.join(consensus)

results = [ consensus , len(consensus) , CS() , SOP() ]

def show():
    print('Column Score Achieved = {:0.3f}'.format(results[2]) )
    print('Sum_Of_Pair Score Achieved = ', results[3])
    print('Length of Consensus Sequence Obtained', results[1])
    print('Consensus Sequence Obtained\n', results[0])

show()
