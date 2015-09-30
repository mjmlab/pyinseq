# pyinseq example data

# example01 - 80 reads across 2 samples

`example01.txt`

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


```
>chr1_0012_04_for
CGACACCGACGACGGTA
>chr1_0012_04_rev
TGTGGTGAAGACCTGTA
>chr1_0012_98_for
AGTTGAACCTCGTCGTA
>chr1_0012_98_rev
TCTCTTCAATGAAGTTA
>chr1_igr0012_for
CAATTTATCTCAACTTA
>chr1_igr0012_rev
CATATATCTAAATTCTA
>chr1_0033_21_for
CAGAACAAGCAGACTTA
>chr1_0033_21_rev
TATAATCTGTTGCTCTA
>chr2_A0492_98_for
ACACCTAAACTCGCATA
>chr2_A0492_98_rev
ACAGGTGCTTCTCGTTA
>chr2_A0492_03_for
ACTGGTTCCAGACATTA
>chr2_A0492_03_rev
AAAAATACATTGGCTTA
>chr2_igrA0492_for
CTATATTACTCACACTA
>chr2_igrA0492_rev
GCAATTTTATGTCAATA
>plas_B0007_96_for
AGATCCAACCTCTTTTA
>plas_B0007_96_rev
AACGGTCTCTTTTGCTA
>plas_B0007_60_for
AATGCAGCACAGATTTA
>plas_B0007_60_rev
CCCAAGAGAAGTGCTTA
>plas_B0007_04_for
GTTGAGGTCTTGGATTA
>plas_B0007_04_rev
TGCACTGGTCGGGAATA
>plas_igrB0007_for
TTCTCTTATGAACATTA
>plas_igrB0007_rev
CACATTAACCTCATTTA
```
