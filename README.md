# Protamine_Exploration
A project in the Cristanti Lab exploring the protamine genes

## UTR Exploration
We first had to extract all UTR sequences around the protamine genes identified using [OrthoExpress](https://github.com/dthorburn/OrthoExpress). However, to not miss any potentially unannotated sequences I created synthetic UTRs by first estiamting the average size of UTRs in this dataset, then adding syntehtic 5' and/or 3' UTRs when the annotation was not present in the GTF files and converts coordiates into BED12 format files: ```Candidate_BED12_SynthUTRs.R```.  
Sequences were then extracted from the BED12 files using the following [BEDTools](https://bedtools.readthedocs.io/en/latest/) command:
```
bedtools getfasta -split -name -s -fi ${ref_fasta} -bed ${BED12} -fo ${Assembly}_Protamine_Candidates_wUTRs.fasta
```
Sequences were then aligned across the dataset using [edialign](https://bioinformaticshome.com/tools/msa/descriptions/edialign.html)
