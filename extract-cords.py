#!/usr/bin/python

import sys
import argparse

def writeArray(fp,name,arr):
    fp.write("const double "+name+'['+str(len(arr))+']={')
    for e in arr:
        if e==0:
            fp.write('0,')
        else:
            fp.write('{0: f},'.format(e))
    fp.write('};\n')

p=argparse.ArgumentParser(description='Extract spherical harmonic coefficients from IGRF file an write to a C file.')
p.add_argument('infile',type=argparse.FileType('r'))
p.add_argument('-o','--outfile',type=argparse.FileType('r'),dest='outfile',default=sys.stdout)

args=p.parse_args()

inf=args.infile
otf=args.outfile


have_header=False

while not have_header:
    l=inf.readline()
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

for line in inf:
    line=line.split()
    g_h=line[0]
    n=line[1]
    m=line[2]

    sv+=[float(line[svIdx])]
    gh+=[float(line[ghIdx])]

inf.close()

writeArray(otf,'gh',gh)
writeArray(otf,'sv',sv)

otf.close()

