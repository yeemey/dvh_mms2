'''
Filter reads that have already been tagged with sample names at fastq header ends.
Sample name tags are preceded by '--'
Usage: python filter_tagged_reads.py filename.fastq.gz R[1 or 2] tag [... tag_n]
'''
from sys import argv
import gzip

filename = argv[1]
r1r2 = '_' + argv[2]
samples = argv[3:]
print(samples)
with gzip.open(filename, 'rt') as file:
    fastq = file.readlines()
    linecount = 0
    while linecount < len(fastq):
        print(fastq[linecount])
        tag = fastq[linecount].split('--')[-1].rstrip()
        if tag in samples:
            read = fastq[linecount:(linecount+4)]
            with open(tag + r1r2 + '.fastq', 'a+') as samplefile:
                for entry in read:
                    samplefile.write(entry)
        linecount += 4
print('Filtered!')
