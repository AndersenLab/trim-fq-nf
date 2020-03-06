#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"

// nextflow main.nf -profile debug
// nextflow main.nf --fastq_folder folder_name 
// The folder_name is the name without full path
// See config in config/ folder for default folder path
// default profile is quest

println "Running fastp trimming on ${params.raw_path}/${params.fastq_folder}"



/* 
    ==================================================
    Calculate MD5 for all files in a single process
    ==================================================
*/

process pre_trim_md5sum {
    // this process runs in the data folder, instead of nextflow working dir

    """
    cd ${params.raw_path}/${params.fastq_folder}
    md5sum *.fq.gz > md5.txt
    """
}

/* 
    ==================
    trim raw data
    ==================
*/

process fastp_trim {

    tag { sampleID }

    publishDir "${params.processed_path}/${params.fastq_folder}", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.processed_path}/${params.fastq_folder}/multi_QC", mode: 'copy', pattern: "*_fastp.html"

    input:
      tuple val(sampleID), path(fq1), path(fq2) 

    output:
      tuple val(sampleID), path("${sampleID}_1P.fq.gz"), path("${sampleID}_2P.fq.gz"), emit: fastq_post
      path "*_fastp.json", emit: fastp_json
      path "*_fastp.html"
    
    """

    fastp -i $fq1 -I $fq2 \\
          -o ${sampleID}_1P.fq.gz -O ${sampleID}_2P.fq.gz \\
          --length_required 20 \\
          -j ${sampleID}_fastp.json -h ${sampleID}_fastp.html

    """
}

/* 
    =======================
    combine all trim report
    =======================
*/

process multi_QC {

    publishDir "${params.processed_path}/${params.fastq_folder}/multi_QC", mode: 'copy'

    input:
      path(json) 

    output:
      path "*.html"
      path "multiqc_data/*"
    
    """

    multiqc . --data-format tsv

    """
}

// read input
fq = Channel.fromFilePairs("${params.raw_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)


workflow { 

    pre_trim_md5sum()  // check sum for all files. 

    fq | fastp_trim

    fastp_trim.out.fastp_json.collect() | multi_QC

}
