# Tutorial

## Requirements

Mac or Linux-based operating systems.
  - See these directions for [running on an Amazon EC2 instance](https://angus.readthedocs.io/en/2016/amazon/index.html).

Python 3.5 or higher (installation described below).

## Install Python

Install the [Anaconda Python 3.5 download](https://www.continuum.io/downloads).  

## Install pyinseq

Install from Github.

```
pip install git+git://github.com/mandel01/pyinseq
```

Test for correct installation

```
make test
```

## Data files needed

**Illumina sequencing reads** for the experiment. File can be uncompressed (.fastq or .fq) or gzip-compressed (.fastq.gz or .fq.gz). [input reads `-i` below]

```
@DGL9ZZQ1:720:C6YD0ACXX:2:1101:1246:2185 1:N:0:
GAAGCGACACCGACGACGGTAACAGGTTGGATGATAAGTCCCCGGTCTTCG
+
CCCFFFDDHHHHHCGHHHIJDHIIJJFHHJIJJJJIJJDHIJJHIAGIJJJ
...
```

**Sample file** describing the sample names and barcodes. Sample names should be restricted to letters, numbers, dash (-), and underscore (_), with a tab between the sample and the barcode in each row of a text file. It is recommended that the file be prepared in a text editor to ensure that additional hidden characters are not introduced. [Textwrangler](http://www.barebones.com/products/TextWrangler/) is recommended for Mac. [sample file `-s` below]

```
E001_01	GAAG
E001_02	CTTT
```

**GenBank file** listing the features and DNA sequence for the organism. If the organism has multiple chromosomes/contigs in the sequence the file they should be concatenated into a single file. Ensure that the double slash `//` at the end of the file remains to separate each contig. [GenBank file `-g` below]

```
LOCUS       CP000020             2897536 bp    DNA     circular BCT 02-APR-2008
DEFINITION  Vibrio fischeri ES114 chromosome I, complete sequence.
ACCESSION   CP000020
VERSION     CP000020.2  GI:171902228
KEYWORDS    .
SOURCE      Vibrio fischeri ES114
...
CDS             complement(313..747)
                /gene="mioC"
                /locus_tag="VF_0001"
...
ORIGIN
        1 aagatcactt aatatatata agatctttta aagagatctt ttattagatc tattatatag
       61 atcgtcgatc tctgtggata agtgataaat gatcaatagg atcatatact ttagatggat
      121 ccaaagttgt tatctttctt tgatcttcga tcggacagct tgaggacaaa agagttagtt
      181 atccacaagg ggggagggcg ttagatctta ttcaatggat aactataact tgatcactgg
      241 atctttctat agttatccac atagtaggta tcatctattt aataactttt atagatcgga
      301 caacactttt tattaacaaa tgtgttgttt tagccacaat tctgctggtt cttcagggat
      361 actattttga gttacatcta tctctaatag agggtacttt tcgatagcgc catgctcttt
      421 taaggtattt tgaagtgatt ttcctgctgc gcagaaagtg tcataacttg aatcaccgat
      481 agcaacaaca gcaaattcaa cgctcgataa tggctgagtc acgctttcta actgttgaat
      541 aaatggttta atgttgtcag ggtaatcacc agcaccgtga gttgatacaa cgagtaacca
      601 taagctatca atatcaatat catctaaatt tggttgatta tgaatatcgg tggaaaaatc
      661 catttcttct aataaatcag caagatggtc accaacatac tcagcaccac ctagagtgct
      721 tcctgtaata atagatactt ttttcatgaa tttatcctat aaaaatataa aaaatgggcc
      781 tacataggcc cattattaat cttattaata ttggttttat ttaccaatac agaatgaagt
...
//
LOCUS       CP000021             1330333 bp    DNA     circular BCT 02-APR-2008
DEFINITION  Vibrio fischeri ES114 chromosome II, complete sequence.
...
```

## pyinseq command line operation

Run the example data set:

```
pyinseq -i <input file> -s <sample file> -g <genbank file> -e <experiment name>
```

`-i` (or `--input`)

- Illumina reads file in FASTQ or gzipped-FASTQ format.

`-s` (or `--samples`)

- Sample list where each line is the `sample_name`tab`barcode(4bp)`.

`-g` (or `--genome`)

- Concatenated Genbank File for all contigs for the genome.

`-e` (or `--experiment`)

- All results files will be created in the subfolder of this name.

## Optional arguments

`-d` (or `--disruption`)

- Five-prime fraction of gene (`0.0` - `1.0`) that must be disrupted for the hit to be counted in the summary_gene_table. Often insertions at the 3' end of a gene do not disrupt function so it may be of interest to run the pipeline with a disruption value of `0.8` or `0.9`.

## Output files

### `results/` directory  

- `summary_gene_table.text`: summary for entire experiment  
- `<sample>_sites.txt` (for each sample): Counts of each insertion in each sample   
- `<sample>_genes.txt` (for each sample): Counts of each insertion mapped to genes  
- `<sample>_bowtie.txt` (for each sample): Bowtie mapping results  
- `<sample>_trimmed.fastq` (for each sample): Demultiplexed fastq reads trimmed for the chromosome sequence only  

### `results/genome_lookup/` subdirectory
- `genome.fna`: genome fasta nucleotide file  
- `genome.ftt`: genome feature file  
- bowtie indexes  

### `results/raw_data/` subdirectory 
- `<sample>.fastq` (for each sample): demultiplexed files for each sample/barcode
- `_other.fastq`: demultiplexed files for unrecognized barcodes
