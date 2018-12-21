#!/usr/bin/python
from Bio import SeqIO
import sys
import os

#usage: python bysize.py in.file out.file format 200 300
format1=sys.argv[3]
input_seq_iterator = SeqIO.parse(open(sys.argv[1], "rU"), format1)
short_seq_iterator = (record for record in input_seq_iterator \
if len(record.seq) > int(sys.argv[4]) and len(record.seq) < int(sys.argv[5]))

output_handle = open(sys.argv[2], "w")
SeqIO.write(short_seq_iterator, output_handle, format1)
output_handle.close()