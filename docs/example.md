# Example data sets

## example01

### 80 reads across 2 samples (plus 2 reads on an untracked sample)

`_exampleData/example01.fastq` In this repository.  

`[_exampleData/example01.txt](../_exampleData/example01.txt)` In this repository.

sampleName  | barcode
------------- | -------------
E001_01  | GAAG
E001_02  | CTTT


feature | index | GAAG reads | GAAG CPM (d 1.0) | GAAG CPM (d 0.9)
---- | ---- | ---- | ----: | ----:
VF_0012 | 0.04 | 2 | 100000 | 50000
VF_0012 | 0.98 | 2 | """""" | .
Intergenic following VF_0012 | -- | 6 | . | .
VF_0033 | 0.21 | 2 | 50000 | 50000
VF_A0492 | 0.03 | 6 | 300000 | 150000
VF_A0492 | 0.98 | 6 | """""" | .
Intergenic following VF_A0492 | -- | 2 | . | .
VF_B0007 | 0.04 | 8 | 300000 | 250000
VF_B0007 | 0.60 | 2 | """""" | """"""
VF_B0007 | 0.96 | 2 | """""" | .
Intergenic following VF_B0007 | -- | 2 | . | .

The CTTT barcode has 2 fewer reads in VF_A0492 (-d 0.03) and 2 additional reads in VF_0033  

The sample set also includes 1 hit (2 reads) that map to VF_0012 on barcode TTTT, which should not be displayed since that barcode is not in the sample list.



## example02

### 100,769,798 reads across 5 samples

`example02.fastq` Download all 100,769,798 reads, [4.09 GB gzipped](http://bit.ly/1MnBq18).  
MD5 checksum (.fastq.gz) = 52a1775768c1c1cb759dd6154b0ec237   
MD5 checksum (.fastq) = 28fa47bc1f6b0b09710e6bd8e8027297  

`example02_AAAA.fastq` Download only the 20,956,625 reads in barcode AAAA, [843 MB gzipped](http://bit.ly/1WXYJWK).  
MD5 checksum (.fastq.gz) = c6975f7f93896250be5302b0af7e8a7b  
MD5 checksum (.fastq) = b7424f831ff6ffc20ba0016968b25c69

`example02.txt` In this repository.

sampleName  | barcode
------------- | -------------
E137_1	| TTTT
E137_2	| AGGA
E137_3	| ATCG
E137_4	| GTCA
E137_5	| AAAA
