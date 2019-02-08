import random

rnd=random.Random()

f=open('myio.py','w')

lk={0:'A',1:'T',2:'G',3:'C'}

for i in range(10):
    f.write(''.join(lk[rnd.randint(0,3)] for i in range(100))+'\n')

f.close()
