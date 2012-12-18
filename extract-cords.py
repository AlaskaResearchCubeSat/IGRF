#!/usr/bin/python

import sys
import argparse
import math

#write a python array to a C source array
def writeArray(fp,vt,name,arr,l=None):
    if l is None:
        l=len(arr)
    #write declaration
    fp.write(vt+" "+name+'['+str(l)+']={')
    #write elements
    for i,e in enumerate(arr[0:l]):
        #check element type
        if type(e) is str:
            #write string element
            fp.write(e)
        else:
            #write double element
            fp.write('{0: .18f}'.format(e))
        #skip comma for last element
        if i!=l-1:
            fp.write(',')
    #write closing bracket 
    fp.write('};\n')

#write variable declaration in C source
def writeVar(fp,vt,name,val):
    fp.write(vt+" "+name+'='+str(val)+';\n')

#initialize parser
p=argparse.ArgumentParser(description='Extract spherical harmonic coefficients from IGRF file an write to a C file.')
p.add_argument('infile',type=argparse.FileType('r'))
p.add_argument('-o','--outfile',type=argparse.FileType('w'),dest='outfile',default=sys.stdout)
p.add_argument('-m','--model-order',type=int,dest='igrf_ord')
p.add_argument('-s','--sv-order',type=int,dest='sv_ord')

#parse arguments
args=p.parse_args()

#flag for header
have_header=False
t=None

#scan lines until all header lines are found
while not have_header:
    #read a line
    l=args.infile.readline().strip()
    #check for comment and skip
    if l.startswith('#'):
        continue
    #split out into columns
    h=l.split()
    #check for header line
    if h[0]=='g/h':
        have_header=True
    #check for type header
    elif h[0] in ['IGRF','DGRF','SV']:
        t=h
    else:
        print("Warning: Unknown line in file {}:\n\t\"{}\"".format(args.infile.name,l))

#find secular variation model in array
i=t.index('SV')
#get length difference in header line and type line for offset
os=len(h)-len(t)
#get index of SV
svIdx=os+i
#get year of SV
svd=int(float(h[svIdx].split('-')[0]))

#magnetic field model parameters are one col before SV
ghIdx=svIdx-1
#get date
date=int(float(h[ghIdx]))


#check that model year matches SV year
if date!=svd:
    print("Error: Secular Variation Date does not match model date")
    exit(1)

#initialize data arrays
gh=[]
sv=[]

for line in args.infile:
    #split into columns
    line=line.split()
    #TODO: make sure that these are correct
    g_h=line[0]
    n=line[1]
    m=line[2]

    #get data from file
    sv+=[line[svIdx]]
    gh+=[line[ghIdx]]

#close input file
args.infile.close()

#Check if base model order given
if args.igrf_ord is None:
    #not given, set defaults
    l=len(gh)
    #find last non zero element
    for i in range(l-1,0,-1):
        if float(gh[i])!=0:
            l=i+1
            break
    #solve order from length
    args.igrf_ord=int(math.sqrt(l+1)-1)
#check if secular variation model order given
if args.sv_ord is None:
    #not given, set defaults
    l=len(sv)
    #find last non zero element
    for i in range(l-1,0,-1):
        if float(sv[i])!=0:
            l=i+1
            break
    #solve order from length
    args.sv_ord=int(math.sqrt(l+1)-1)

#includes
args.outfile.write("#include \"igrfCoeffs.h\"\n\n")
#write date
writeVar(args.outfile,'const int','igrf_date',date);
#write sv order
writeVar(args.outfile,'const int','igrf_ord',args.igrf_ord);
#write model  order
writeVar(args.outfile,'const int','sv_ord',args.sv_ord);
#write coefficients
writeArray(args.outfile,'const float','igrf_coeffs',gh,args.igrf_ord*(args.igrf_ord+2))
#add some new lines between
args.outfile.write("\n\n")
#write specular variations
writeArray(args.outfile,'const float','igrf_sv',sv,args.sv_ord*(args.sv_ord+2))
#add some lines at the end
args.outfile.write("\n\n")

#print completion message if not writing to stdout
if args.outfile is not sys.stdout:
    #close outfile
    args.outfile.close()
    #print completion message
    print("File \"{}\" written successfully".format(args.outfile.name))

