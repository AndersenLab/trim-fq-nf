#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
*/
nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"

include md5sum as md5sum_pre from './md5.module.nf'
include md5sum as md5sum_post from './md5.module.nf'

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

md5sum_path = "${params.processed_path}/${params.fastq_folder}/md5sums.txt"

workflow { 
      
    fq | (md5sum_pre & fastp_trim)
    fastp_trim.out.fastq_post | md5sum_post // Run md5sum on post-trim
    md5sum_pre.out.concat(md5sum_post.out).collectFile(name: md5sum_path, newLine:false)

    fastp_trim.out.fastp_json.collect() | multi_QC

}
