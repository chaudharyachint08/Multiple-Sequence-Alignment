import numpy as np

dct  = { '-': 0  , 'A': 1  , 'T': 2  , 'G': 3  ,  'C': 4  }
dct2 = {  0 :'-' ,  1 :'A' ,  2 :'T' ,  3 :'G' ,   4 :'C' }

file = open(input('Enter File Containing DNA Sequences to be aligned using AC0-S\n'))

NS,LN = [int(x) for x in input('Enter Number of DNA Sequences & length to be taken into consideration\n').split()]

mat = [[j for j in file.readline().strip()[:LN]]for i in range(NS)]


def lcs1(X , Y):
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

NC,AC,SL,EL,DD = [int(x) for x  in input('Enter No. of Cycles, Pop_per_Cycle,Start_Len,End_Length, Drift\n').split()]
RND, EVAP, INTENSITY = [float(x) for x  in input('Enter Probability of Random_Chance, Rate of Evaporation & Intensity of Trail\n').split()]


class Itrail:
    '''ITrail class to maintain Trail associated with each sequence'''

    
    def __init__(self,const,seq):
        const = float(const)
        self.trail = np.array([const]*LN)
        self.seq   = seq

    def evap(self):
        self.trail*= (1-EVAP/100)


class Ant:
    '''Each Ant is associated with a Itrail'''

    def __init__(self,trail,length,pos):
        self.length = length
        self.pos    = pos
        self.seq    = trail.seq
        self.trail  = trail

    def make_drift(self,other):
        all_pos = [ x for x in range(self.pos-DD , self.pos+DD+1) if x>=0 and x<(LN-self.length+1) ]
        k,val = -1,-1
        for i in all_pos:
            new_val = lcs1( self.seq[self.pos:self.pos+self.length] , other.seq[i:i+self.length])
            if new_val > val:
                k = i
                val = new_val

        if np.random.rand() < RND:
            k = np.random.choice(all_pos)
            
        other.trail[ k : k+self.length ] += INTENSITY


combined = [Itrail(0,mat[i]) for i in range(NS)]

rate_of_length_increase = (EL-SL)/(NC-1)

for i in range(NC):
    length = int(round(SL + rate_of_length_increase*i))
    valid_pos = [x for x in range(LN-length+1)]
    
    Ant_Colony = [[Ant(combined[i],length,np.random.choice(valid_pos)) for j in range(AC)] for i in range(NS)]

    for i in range(NS):
        for j in range(NS):
            if i!=j:
                for k in range(AC):
                    Ant_Colony[i][k].make_drift(combined[j])

    for i in range(NS):
        combined[i].evap()

#for i in range(NS):
#    print(combined[i].trail)

sequences_with_gaps = [''.join( ( (i.seq[k] if i.trail[k]>=0.5*i.trail.max() else '-') for k in range(LN) ) ) for i in combined]

consensus = [[(k,sum(combined[i].trail[j] for i in range(NS) if combined[i].seq[j]==k)) for k in ('A','T','G','C')] for j in range(LN)]
    
final_sequence = []
for i in consensus:
    i.sort(key=lambda x:x[1],reverse=True)
    if i[0][1] > i[1][1]*1.5:
        final_sequence.append(i[0][0])
    
results = [ ''.join(final_sequence) ]


def CS():

    score = 0

    for i in range(LN):
        su = { 'A':0 , 'T':0 , 'G':0 , 'C':0 , '-':0 }

        for j in sequences_with_gaps:
            su[j[i]] += 1

        m = max(su.values())
        score += (m*(1+m/LN))
    
    return score


def SOP():
    '''Sum_Of_Pair Function'''

    score = 0

    nA,nB,nC = 0,0,0
    A,B,C    = 9,-5,-2

    for i in sequences_with_gaps:
        for j in sequences_with_gaps:
            if i is not j:
                for k in range(LN):
                    if i[k]==j[k]:
                        nA += 1
                    elif '-' in (i[k],j[k]):
                        nB += 1
                    else:
                        nC += 1

    score = A*nA - B*nB - C*nC

    return str(nA)+'*A - '+str(nB)+'*B - '+str(nC)+'*C' , score


results.extend([len(results[0]),CS(),SOP()])

def show():
    print('Column Score Achieved = {:0.3f}'.format(results[2]) )
    print('Sum_Of_Pair Score Achieved = ', results[3])
    print('Length of Consensus Sequence Obtained', results[1])
    print('Consensus Sequence Obtained\n', results[0])

show()
