#!/usr/bin/python

from Bio import SeqIO
import sys
import os


#Theo Allnutt 2016
#Converts files of different sequence format
#usage: convert.py infile informat outformat [see http://biopython.org/wiki/SeqIO] 

inputfile = sys.argv[1]
informat = str(sys.argv[2])
outformat = str(sys.argv[3])
try:
	num = str(sys.argv[4])
except:
	num=0
	
f = open(inputfile,'rb')
outputfile = inputfile[:len(inputfile)-len(informat)] + outformat
g = open(outputfile,'w')
if num<>0:
	count = SeqIO.convert(f[0:num], informat,g, outformat)
else:
	count = SeqIO.convert(f, informat,g, outformat)

print("Converted %i records" % count)
