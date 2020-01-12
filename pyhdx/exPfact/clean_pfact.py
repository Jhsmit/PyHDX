import os

list1=[x.strip() for x in open('none.dat', 'r').readlines()]

for f in os.listdir():
    if f.endswith('.pfact'):
        lines = open(f, 'r').readlines()
        fout = open(f.split('.p')[0]+'.cpfact', 'w')
        for line in lines:
            if not line.strip().split()[0] in list1:
                fout.write(line)
        fout.close()
