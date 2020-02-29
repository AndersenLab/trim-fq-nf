#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Dan Lu <dan.lu@northwestern.edu>
*/


nextflow.preview.dsl=2


params.fastq_path="/projects/b1059/data/fastq/WI/dna/raw"
params.trimmed_path="/projects/b1059/data/fastq/WI/dna/processed"

params.fastq_folder="test_data/fastq"
params.trimmed_folder="${params.fastq_folder}_fastp"




/* 
    ==============================================
    FASTQC on raw data and combine with MultiQC
    ==============================================
*/


process pre_trim_fastqc {

    tag { sampleID }

    publishDir "${params.fastq_path}/${params.fastq_folder}/fastqc", mode: 'copy', pattern: "*.html"

    input:
      tuple sampleID, path(fq1), path(fq2) 

    output:
    /* note each item in the channel is a file pair */
       path "*.html" // output to fastqc folder
       path "*.zip", emit: pre_trim_fastqc_zip // emit for the next process

    """
    fastqc $fq1
    fastqc $fq2
    """
}




process pre_trim_multi_QC {


    publishDir "${params.fastq_path}/${params.fastq_folder}/fastqc", mode: 'copy', pattern: "*.html"

    input:
      path(fastqc_zip)

    output:
      path "*.html"
    
      """

      multiqc .

      """
}


/* 
    ==================
    trim raw data
    ==================
*/


process fastp_trim {


    publishDir "${params.trimmed_path}/${params.trimmed_folder}", mode: 'copy', pattern: "*.fq.gz"


    input:
      tuple sampleID, path(fq1), path(fq2) 

    output:
      path "*_trimmed.fq.gz" 
      path "*_fastp.json", emit: fastp_json
    

      """

      fastp -i $fq1 -I $fq2 \\
            -o ${sampleID}_1_trimmed.fq.gz -O ${sampleID}_2_trimmed.fq.gz \\
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


    publishDir "${params.trimmed_path}/${params.trimmed_folder}/multi_QC", mode: 'copy', pattern: "*.html"

    input:
      path(json) 

    output:
      path "*.html"
    
    script:
      """

      multiqc .

      """
}




// read input
fq = Channel.fromFilePairs("${params.fastq_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)


// run workflow
workflow { 

fq | (pre_trim_fastqc & fastp_trim)

pre_trim_fastqc.out.pre_trim_fastqc_zip.flatten().toSortedList() | pre_trim_multi_QC

fastp_trim.out.fastp_json.toSortedList() | multi_QC

}