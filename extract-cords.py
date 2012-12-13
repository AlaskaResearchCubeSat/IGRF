#!/usr/bin/python

import sys

def writeArray(fp,name,arr):
    fp.write("const double "+name+'['+str(len(arr))+']={')
    for e in arr:
        if e==0:
            fp.write('0,')
        else:
            fp.write('{0: f},'.format(e))
    fp.write('};\n')


fname='igrf11coeffs.txt'

fp=open(fname,'r')

have_header=False

while not have_header:
    l=fp.readline()
    if l.startswith('#'):
        continue
    h=l.split()
    if h[0]=='g/h':
        have_header=True

#initialize coordinate array

svIdx=26
ghIdx=svIdx-1
gh=[]
sv=[]

for line in fp:
    line=line.split()
    g_h=line[0]
    n=line[1]
    m=line[2]

    sv+=[float(line[svIdx])]
    gh+=[float(line[ghIdx])]

fp.close()

writeArray(sys.stdout,'gh',gh)
writeArray(sys.stdout,'sv',sv)

