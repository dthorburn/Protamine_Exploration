# Protamine_Exploration
A project in the Cristanti Lab exploring the protamine genes

## UTR Exploration
We first had to extract all UTR sequences around the protamine genes identified using [OrthoExpress](https://github.com/dthorburn/OrthoExpress). However, to not miss any potentially unannotated sequences I created synthetic UTRs by first estiamting the average size of UTRs in this dataset, then adding syntehtic 5' and/or 3' UTRs when the annotation was not present in the GTF files and converts coordiates into BED12 format files: ```Candidate_BED12_SynthUTRs.R```.  
Sequences were then extracted from the BED12 files using the following [BEDTools](https://bedtools.readthedocs.io/en/latest/) command:
```
bedtools getfasta -split -name -s -fi ${ref_fasta} -bed ${BED12} -fo ${Assembly}_Protamine_Candidates_wUTRs.fasta
```

## Phylogenetic Tree Construction
To place everything into an evolutionary context, we then aligned the sets of sequences with [MUSCLE](https://www.ebi.ac.uk/Tools/msa/muscle/), [MAFFT](https://mafft.cbrc.jp/alignment/software/), and [edialign](https://bioinformaticshome.com/tools/msa/descriptions/edialign.html). 

I've developed a [nextflow](https://www.nextflow.io/) pipeline (version 20.10.0; Imperial HPC's most recent version 03/02/2022) to take multiple sequence fasta files, align them with the choice of either MUSCLE or MAFFT (edalign was omitted due to module availablity, but the pipeline can be updated to include an edalign conda environmente easily). This pipeline was developed to use the modules available on Imperial College's HPC. If you don't have access you'll need to update how the software is loaded into the environment. Phylogenetic trees are constructed using RAxML. Use the `--version` option to see the versions used during development. Aligned sequences can be provided and used with the `--Skip_Aln` option. 

### Usage
To use this pipeline follow these instructions:

1. Download the `RAxML_Tree.*` files into your project directory. 
2. Update the `RAxML_Tree.config` file to reflect your needs. For the majority of cases, you'll only need to update the input directory path (i.e., `UnalignDir` if you have not aligned the sequences, or `AlignDir` if you are providing aligned sequences). *NB. The pipeline will create a tree for each fasta file in the directory*. 
3. Update the project directory path in the `RAxML_Tree.sh` file. 
4. Run the pipeline using the command `qsub RAxML_Tree.sh`

Below is the help message from `RAxML_Tree.nf`:
```
Usage:
  This pipeline was developed to use on Imperial College London's HPC. Alternatively, you can set up a conda environment with
  all the tools, but you'd have to alter the process beforeScripts to set up the environment appropriately.

  First, alter the RAxML_Tree.config file to reflect the paths you'll be using. Best practise is to leave the paths of all but
  the inital path to your input files. The pipeline expects fasta files that can be detected with the glob *.fa{,sta,a}. i.e.,
  .fasta, .faa, .fa.

  Once the config file is updated, run with the following command run from the project directory:
  qsub RAxML.sh

  Directory Structure:
    /Project_dir/
      | - RAxML.config                                      Config file to update - Required
      | - RAxML.nf                                  
      | - RAxML.sh                                          Nextflow job submitter. Update if skipping steps or adding paramaters
      | - 01_Aligned/
      | - 02_Tree_Parts/
      | - 03_Final_Trees/                                   Final results will be located here. 
  
  Optional Paramaters:
    --help                                                  Show this message
    --version                                               See versions used to develop pipeline
    --Skip_Aln                                              Skips the fasta alignment step
    --Skip_ML                                               Skips creating the ML tree
    --Skip_BS                                               Skips creating the bootstraps
    --Skip_BuildTree                                        Skips the final tree construction
    --Aln_Mode "muscle"                                     Choice of alignment software (muscle/mafft; default muscle)    
```
