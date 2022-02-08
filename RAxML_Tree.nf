#!/usr/bin/env nextflow

// Pipeline developed for building a ML phylogenetic tree based upon a fasta file. 
// Author: Miles Thorburn <d.thorburn@imperial.ac.uk>
// Date last modified: 02/02/2022

def helpMessage() {
  log.info """
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
           """
}
def versionMessage() {
  log.info  """
            For repeatability, here are the versions I used to construct the pipeline:
              muscle/3.7
              mafft/7.271
              raxml/8.2.9
            """  
}

println "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nPhylogenetic Tree Pipe v0.1\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
// println "${PWD}, ${HOME}, ${PATH}"

// Show help message
if (params.help) {
    helpMessage()
    exit 0
}

// Show version information
if (params.version) {
    versionMessage()
    exit 0
}

// Step 1: ALIGNMENT
if( params.Skip_Aln == false ) {
  Channel
     .fromPath("${params.UnalignDir}/*.fa{,sta,a}")
     .ifEmpty { error "Cannot find input folder ${params.UnalignDir}" }
     .set { raw_fasta_ch }

  process Alignment {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.Aln_Forks

    publishDir(
      path: "${params.AlignDir}",
      mode: 'copy',
    )

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.ALN_threads}:mem=${params.ALN_memory}gb -lwalltime=${params.ALN_walltime}:00:00"

    input:
    path fasta from raw_fasta_ch

    output:
    tuple val(SampleID), path("${SampleID}_aligned.fasta") into ml_in_ch, bs_in_ch


    beforeScript 'module load mafft/7.271; module load muscle/3.7'
    
    script:
    // Getting the file basename
    SampleID = fasta.baseName

    if( params.Aln_Mode == 'mafft' )
      """
      mafft --thread ${params.ALN_threads} ${params.ALN_mafft_argmts} --maxiterate ${params.ALN_mafft_Iters} ${fasta} > ${SampleID}_aligned.fasta
      """
    else if( params.Aln_Mode == 'muscle' )
      """
      ALN_muscle_mb=`expr ${params.ALN_memory} \\* 1000`
      muscle -maxmb \${ALN_muscle_mb} -in ${fasta} -out ${SampleID}_aligned.fasta
      """
    else
      error "Invalid alignment mode: ${params.Aln_Mode}, use either mafft or muslce"
  }
}

// Step 2: ML Tree
if( params.Skip_ML == false ) {
  // Handles cases where alignment is skipped. Loads the tuple files into the right channel
  if( params.Skip_Aln ){
    Channel
      .fromPath("${params.AlignDir}/*.fa{,sta,a}")
      .ifEmpty { error "Cannot find input folder ${params.AlignDir}" }
      .map { file -> tuple(file.baseName, file) }
      .set { ml_in_ch }
  }

  process ML_Tree {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.ML_Forks

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.ML_threads}:mem=${params.ML_memory}gb -lwalltime=${params.ML_walltime}:00:00"

    publishDir(
      path: "${params.TreePartsDir}",
      mode: 'copy',
      saveAs: { filename -> "${SampleID}.ML" }
    )

    input:
    set SampleID, path(aln_fas) from ml_in_ch

    output:
    tuple val (SampleID), path("RAxML_bestTree.${SampleID}.nwk") into ml_tree

    beforeScript 'module load raxml/8.2.9'

    script:
    """
    raxmlHPC-AVX -m ${params.ML_Model} -p ${params.ML_Seed} -N ${params.ML_Iters} -s ${aln_fas} ${params.ML_argmts} -n ${SampleID}.nwk
    """
  }
}

// Step 3: Bootstrapping Tree
if( params.Skip_BS == false ){
  // Handles cases where alignment was skipped by loading tuple files into the correct channel
  if( params.Skip_Aln ){
    Channel
      .fromPath("${params.AlignDir}/*.fa{,sta,a}")
      .ifEmpty { error "Cannot find input folder ${params.AlignDir}" }
      .map { file -> tuple(file.baseName, file) }
      .set { bs_in_ch }
  }

  // Generating the channel to parallelise this step
  Nums = Channel.of(1..params.BT_par_jobs)
  // Combining the input channel. Important which comes first in how the value is added to the tuple, here it is "file.baseName, file, Num"
  // .combine is needed to iterate input channel, otherwise it will simply merge the tuples
  bs_in_ch
    .combine(Nums)
    .set { Parallel_bs_in_ch }
   
  process Bootstrapping {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 600 as long); return 'retry' }
    maxRetries 3
    maxForks params.Boot_Forks

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BS_threads}:mem=${params.BS_memory}gb -lwalltime=${params.BS_walltime}:00:00"

    // Deprecated due to output channel being consumed by collectFile operator which publishes the file
    //publishDir(
    //  path: "${params.TreePartsDir}",
    //  mode: 'copy',
    //)

    input:
    set SampleID, path(aln_fas), rep from Parallel_bs_in_ch

    output:
    tuple val(SampleID), path("RAxML_bootstrap.${SampleID}_bootstraps.${rep}") into raw_bs_trees

    beforeScript 'module load raxml/8.2.9'

    // bash doesn't emit floating points, so incorrect division will result in a decrease in total bootsraps. 
    // Added a check to ensure the parallelisation will result in the correct number of bootstraps being generated. 
    script:
    """
    seed=\$RANDOM
    echo \$seed 
    reps_here=`expr ${params.BS_bootstraps} / ${params.BT_par_jobs}`
    just_checking=`expr \${reps_here} \\* ${params.BT_par_jobs}`
    
    if [ ${params.BS_bootstraps} == \${just_checking} ]
    then
      raxmlHPC-AVX -f d -m ${params.ML_Model} -b \$seed -p \$seed -N \$reps_here ${params.BS_argmts} -s ${aln_fas} -n "${SampleID}_bootstraps.${rep}"
    elif [ ${params.BS_bootstraps} != \${just_checking} ]
    then
      echo "ERROR: Incorrect parallelisation, only \${just_checking} bootsraps will be generated and not ${params.BS_bootstraps}"
      exit 1
    fi
    """
  }
  // Collects all the parallel bootstrap files and puts them into a single file for each SampleID. Uses the first value in the tuple as the grouping criteria. 
  // This should make it so this process runs afer all the bootstrapping is complete. 
  raw_bs_trees
    .collectFile(storeDir: "${params.TreePartsDir}") { SampleID, trees -> 
      [ "${SampleID}.BOOTS",  trees ]}
    .set { all_trees }
}

// Step 4: Building Bipartitions
if( params.Skip_BuildTree == false ){
  // Handling when either the ML or BS steps are skipped. Mostly for debugging. 
  if( params.Skip_ML ){
    Channel
      .fromPath("${params.TreePartsDir}/*.ML")
      .ifEmpty { error "No ML tree detected in $ProjectDir/02_Tree_Parts" }
      .map { file -> tuple(file.baseName, file) }
      .set { ml_tree }
  }
  if( params.Skip_BS ){
    //Channel
    //  .fromPath("${params.TreePartsDir}/RAxML_bootstrap.*")
    //  .ifEmpty { error "No bootstraps detected in $ProjectDir/02_Tree_Parts" }
    //  .map { file -> tuple(file.baseName, file) }
    //  .set { raw_bs_trees }
    Channel
        .fromPath("${params.TreePartsDir}/*.BOOTS")
        .ifEmpty { error "No bootstrap file detected in ${params.TreePartsDir}" }
        .map { file -> tuple(file.baseName, file) }
        .set { all_trees }
  }

  // Combining the files into a single tuple to ensure the correct files are used for each sample
  // Using the 'by' operator means the first object in the tuple of each channel is used as the index. 
  // Indices are zero-based. Do I need to add .collect()?
  ml_tree
    .combine( all_trees, by: 0 )
    .set { all_tree_parts }

  process BuildTree {
    errorStrategy { sleep(Math.pow(2, task.attempt) * 1800 as long); return 'retry' }
    maxRetries 3
    maxForks params.Build_Forks

    executor = 'pbspro'
    clusterOptions = "-lselect=1:ncpus=${params.BT_threads}:mem=${params.BT_memory}gb -lwalltime=${params.BT_walltime}:00:00"

    publishDir(
      path: "${params.OutDir}",
      mode: 'move',
    )

    beforeScript 'module load raxml/8.2.9'

    input:
    //set SampleID_ml, path(ml_tree) from ml_tree
    //set SampleID_bs, path(bs_trees) from all_trees
    set SampleID, path(ml_tree), path(bs_trees) from all_tree_parts

    output:
    path("*")

    script:
    """
    raxmlHPC-AVX -p ${params.BT_Seed} -m GTRCAT -f b ${params.BT_argmts} -t ${ml_tree} -z ${bs_trees} -n ${SampleID}_Final_Supported.nwk
    """
  }
}

