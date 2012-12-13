#!/usr/bin/python

import sys
import argparse

def writeArray(fp,vt,name,arr):
    fp.write(vt+" "+name+'['+str(len(arr))+']={')
    for e in arr:
        if type(e) is str:
            fp.write(e+',')
        else:
            if e==0:
                fp.write('0,')
            else:
                fp.write('{0: .18f},'.format(e))
    fp.write('};\n')

def writeVar(fp,vt,name,val):
    fp.write(vt+" "+name+'='+str(val)+';\n')

p=argparse.ArgumentParser(description='Extract spherical harmonic coefficients from IGRF file an write to a C file.')
p.add_argument('infile',type=argparse.FileType('r'))
p.add_argument('-o','--outfile',type=argparse.FileType('w'),dest='outfile',default=sys.stdout)

args=p.parse_args()

inf=args.infile
otf=args.outfile

date=2010

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
gh=[0]
sv=[0]

for line in inf:
    line=line.split()
    g_h=line[0]
    n=line[1]
    m=line[2]

    sv+=[line[svIdx]]
    gh+=[line[ghIdx]]

inf.close()


otf.write("#include \"igrfCoeffs.h\"\n\n")

#write date
writeVar(otf,'const int','igrf_date',date);
#write sv order
writeVar(otf,'const int','igrf_ord',13);
#write model  order
writeVar(otf,'const int','sv_ord',8);
#write coefficients
writeArray(otf,'const double','igrf_coeffs',gh)
#add some new lines between
otf.write("\n\n")
#write specular variations
writeArray(otf,'const double','igrf_sv',sv)
#add some lines at the end
otf.write("\n\n")

otf.close()

