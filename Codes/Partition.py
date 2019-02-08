f1 = open('DNA.txt')
f2 = open('IP.py','w')

st = f1.readline()

f2.write('X='+'\''+ st[            : len(st)//2 ] +'\''+'\n')
f2.write('Y='+'\''+ st[ len(st)//2 : len(st)    ] +'\''+'\n')

f1.close()
f2.close()
