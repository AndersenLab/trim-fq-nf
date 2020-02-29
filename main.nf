#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Dan Lu <dan.lu@northwestern.edu>
*/

params.fastq_path="/projects/b1059/data/fastq/WI/dna/raw"
params.trimmed_path="/projects/b1059/data/fastq/WI/dna/processed"

params.fastq_folder="test_data/fastq"
params.trimmed_folder="${params.fastq_folder}_fastp"


fq=Channel.fromFilePairs("${params.fastq_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)



process fastp_trim {

    conda "fastp"

    publishDir "${params.trimmed_path}/${params.trimmed_folder}", mode: 'copy', pattern: "*.fq.gz"
/*    publishDir "${params.trimmed_path}/${params.trimmed_folder}", mode: 'copy', pattern: "*_fastp.json"
    publishDir "${params.trimmed_path}/${params.trimmed_folder}", mode: 'copy', pattern: "*_fastp.html"
*/

    input:
      set val(sampleID), file(fq1), file(fq2) from fq

    output:
      set file("*_trimmed.fq.gz"), file("*_fastp.json"), file("*_fastp.html") into fq_trimmed
      file("*_fastp.json") into trim_json
    
    script:
      """

      fastp -i $fq1 -I $fq2 \\
            -o ${sampleID}_1_trimmed.fq.gz -O ${sampleID}_2_trimmed.fq.gz \\
            --length_required 20 \\
            -j ${sampleID}_fastp.json -h ${sampleID}_fastp.html

      """
}



process multi_QC {

    conda "multiqc"

    publishDir "${workflow.launchDir}/multi_QC", mode: 'copy', pattern: "*.html"

    input:
      file(json) from trim_json.toSortedList()

    output:
      file("*.html") into multiQC_report
    
    script:
      """

      multiqc .

      """
}
