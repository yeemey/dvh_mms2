'''
Filter reads that have already been tagged with sample names at header end.
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
