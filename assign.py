#!/usr/bin/env python
"""
Assigns barcodes and transposon data, trims sequence and quality

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.
Unrecognized barcodes get written to file_Other.fastq

"""

import gzip
import sys
import os
import argparse

def parseArgs(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',
        help='input Illumina reads file',
        required=True)
    parser.add_argument('-s', '--samples',
        help='sample list with barcodes',
        required=True)
    parser.add_argument('-e', '--experiment',
        help='experiment name (no spaces or special characters)',
        required=True)
    return parser.parse_args(args)

def assignAndTrim(fastqFile, sampleFile, experiment, tempDir):
    if fastqFile.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    # barcodes = list of barcodes (sequences only)
    # bLen = barcode length

    with open(sampleFile, 'r') as inputSampleFile:
        barcodes, bLen = barcodesPrep(inputSampleFile)

    # initialize to count reads per barcode
    countList = {}
    for b in barcodes:
        countList[b] = 0

    with opener(fastqFile, 'r') as f:

        # truncates the file with the 'w' option before writing data below
        # in cases where the output file already exists
        with open('{0}{1}_Assigned.fastq'.format(tempDir, experiment), 'w') as fo:

            print('Assigning barcode and transposon information.')

            # fastq record
            # identifier =     record[0]
            # sequence =       record[1]
            # annotations =    record[2]
            # quality =        record[3]
            record = []
            for i,line in enumerate(f):
                record.append(line.rstrip('\n'))
                if i % 4 == 3:
                    # barcode sequence in read
                    bc = record[1][0:bLen]

                    # Tn location as number of nuceotides after the barcode
                    tnLoc = record[1].find('ACAGGTTG') - bLen # Tn location

                    # (L) Left / (R) Right / (I) Identical
                    # No code yet for L/R
                    tnEnd = 'I'

                    # Insertion at a TA dinucleotide
                    ta = 'N'
                    if (record[1][(tnLoc+bLen-2):(tnLoc+bLen)]) == "TA":
                        ta = 'Y'

                    # Write in FASTQ format
                    # @ID//experiment:barcode:tnLoc_ta(Y/N)_tnEnd(I/L/R)

                    # Quality filter:
                    # 16-17 bp of putative chromosomal DNA with terminal TA
                    if ta == 'Y' and tnLoc in range(16, 18):

                        # Write only if passes filtering
                        fo.write('{0}//{1}:{2}:{3}_{4}_{5}\n{6}\n{7}\n{8}\n'.format(record[0],
                            experiment,bc,tnLoc,ta,tnEnd,
                            record[1][bLen:(tnLoc+bLen)],
                            record[2],
                            record[3][bLen:(tnLoc+bLen)]
                            ))

                        # Count only if passes filtering
                        if bc in barcodes:
                            countList[bc] += 1

                    record = []

                if ((i+1.0)/4) % 1E6 == 0:
                    print('Reads processed: {:,}'.format(int((i+1)/4)))

            print('Reads processed: {:,}'.format(int((i+1)/4)))

            filteredSum = 0
            for b in countList:
                print('{0}\t{1:,}'.format(b,countList[b]))
                filteredSum = filteredSum + countList[b]
            print('{0:,} ({1:.2%}) reads passed filtering for length and TA terminus'.format(filteredSum, float(filteredSum)/((i+1)/4)))

        # Return output file name for bowtie call

# ===== Start here ===== #

def main():
    #args = parseArgs(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    #assignAndTrim(args.input, args.samples)
    experiment = 'E1'

    createDirectories(experiment)

if __name__ == '__main__':
    main()
