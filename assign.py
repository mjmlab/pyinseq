#!/usr/bin/env python
"""
Assigns barcodes and transposon data, trims sequence and quality

Output file includes the input file name and the barcode:
e.g., file.fastq demultiplexes into file_CGAT.fastq, etc.
Unrecognized barcodes get written to file_Other.fastq

"""

import sys
import os
import argparse

def parse_args(args):
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

def create_directories(experiment):
    """
    Create the project directory and subdirectories

    Attempt to create the directory structure:
    /[experiment]
        /genome
        /temp

    If /experiment directory already exists exit and return error message and
    the full path of the present directory to the user"""

    # Check that experiment name has no special characters or spaces
    pass

    # ERROR MESSAGES
    _error_directory_exists = \
    'PyINSeq Error: The directory already exists for experiment {0}\n' \
    'Delete or rename the {0} directory, or provide a new experiment\n' \
    'name for the current analysis'.format(experiment)

    # Create path or exit with error if it exists.
    try:
        os.makedirs('{}/temp/'.format(experiment))
    except OSError:
        print(_error_directory_exists)
        exit(1)


def barcodes_prep(samples):
    """
    Extract barcodes from an INSeq sample list and conduct basic quality checks

    Input samples is a tab-delimited file. On each line the sample identifier
    should be listed (no spaces), then after a tab the DNA barcode.
    Comment lines begin with '#' and are not read. E.g.:
    # Experiment E01 sample file by M. Mandel
    Input1<tab>CGAT
    Input2<tab>GCTA
    The function checks that all barcodes are the same length.
    The function returns a list of barcodes and sample names and the
    length of the barcodes.

    """

    # Extract barcode list from tab-delimited sample file.
    # Ensure uppercase and stripped
    print('Checking barcodes:')
    barcode_list = {}
    for line in samples:
        if not line.startswith('#'):
            new_barcode = line.rstrip().split('\t')[1]
            if new_barcode in barcode_list:
                print('Error: redundant barcode {}'.format(new_barcode))
                exit(1)
            new_sample = line.rstrip().split('\t')[0]
            if new_sample in barcode_list.values():
                print('Error: redundant sample identifier {}'.format(new_sample))
                exit(1)
            barcode_list[new_barcode] = new_sample
            barcode_length = len(new_barcode)

    # Print barcodes and check for same length
    for b in sorted(barcode_list):
        print('{0}\t{1}'.format(b, barcode_list[b]))
        if len(b) != barcode_length:
            print('Error: barcodes are not the same length')
            exit(1)

    print('n={0} unique barcodes of length {1} nt'.format(len(barcode_list), barcode_length))

    return barcode_list.keys(), barcode_length

def assign_and_trim(fastq_file, sample_file, experiment, temp_dir):
    if fastq_file.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    # barcodes = list of barcodes (sequences only)
    # b_len = barcode length

    with open(sample_file, 'r') as input_sample_file:
        barcodes, b_len = barcodes_prep(input_sample_file)

    # initialize to count reads per barcode
    count_list = {}
    for b in barcodes:
        count_list[b] = 0

    with opener(fastq_file, 'r') as f:

        # Not needed anymore now that the experiment name gets passed in
        """experiment = 'Exp001'   # supply in command line arguments; check alphanumeric only
        if experiment:
            if not experiment.isalnum():
                print('Error: Experiment name should be alphanumeric only. You entered {}'.format(file_prefix))
                exit(1)
            file_prefix = experiment + '_'
            print(file_prefix)"""

        # truncates the file with the 'w' option before writing data below
        # in cases where the output file already exists
        with open('{0}{1}_assigned.fastq'.format(temp_dir, experiment), 'w') as fo:

            print('Assigning barcode and transposon information.')

            # fastq record
            # identifier =     record[0]
            # sequence =       record[1]
            # alt_identifier = record[2]
            # quality =        record[3]
            record = []
            for i,line in enumerate(f):
                record.append(line.rstrip('\n'))
                if i % 4 == 3:
                    # barcode sequence in read
                    bc = record[1][0:b_len]

                    # Tn location as number of nuceotides after the barcode
                    tn_loc = record[1].find('ACAGGTTG') - b_len # Tn location

                    # (L) Left / (R) Right / (I) Identical
                    # No code yet for L/R
                    tn_end = 'I'

                    # Insertion at a TA dinucleotide
                    ta = 'N'
                    if (record[1][(tn_loc+b_len-2):(tn_loc+b_len)]) == "TA":
                        ta = 'Y'

                    # Write in FASTQ format
                    # @ID//experiment:barcode:tnloc_ta(Y/N)_tnend(I/L/R)

                    # Quality filter:
                    # 16-17 bp of putative chromosomal DNA with terminal TA
                    if ta == 'Y' and tn_loc in range(16, 18):

                        # Write only if passes filtering
                        fo.write('{0}//{1}:{2}:{3}_{4}_{5}\n{6}\n{7}\n{8}\n'.format(record[0],
                            experiment,bc,tn_loc,ta,tn_end,
                            record[1][b_len:(tn_loc+b_len)],
                            record[2],
                            record[3][b_len:(tn_loc+b_len)]
                            ))

                        # Count only if passes filtering
                        if bc in barcodes:
                            count_list[bc] += 1

                    record = []

                if ((i+1.0)/4) % 1E6 == 0:
                    print('Reads processed: {:,}'.format(int((i+1)/4)))

            print('Reads processed: {:,}'.format(int((i+1)/4)))

            filtered_sum = 0
            for b in count_list:
                print('{0}\t{1:,}'.format(b,count_list[b]))
                filtered_sum = filtered_sum + count_list[b]
            print('{0:,} ({1:.2%}) reads passed filtering for length and TA terminus'.format(filtered_sum, float(filtered_sum)/((i+1)/4)))

        # Return output file name for bowtie call

# ===== Start here ===== #

def main():
    #args = parse_args(sys.argv[1:])
    #print(args.input)
    #print(args.samples)
    #assign_and_trim(args.input, args.samples)
    experiment = 'E1'

    create_directories(experiment)

if __name__ == '__main__':
    main()
