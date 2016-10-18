## Makefile to run the example script and the test script

## -i / --input, Illumina reads file in FASTQ format
INPUT = data/example/example01.fastq

## -s / --samples, sample list where each line is the sample_name-tab-barcode(4bp)
SAMPLES = data/example/example01.txt

## -g / --genome, concatenated genbank file for all contigs for the genome
GENBANK = data/example/ES114v2.gb

## -e / --experiment, all results files will be created in the subfolder of this name
EXPERIMENT = example01

## Location of test script
TEST = -m pytest pyinseq/tests/

default:
	python run.py -i $(INPUT) -s $(SAMPLES) -g $(GENBANK) -e $(EXPERIMENT)
	python $(TEST)
