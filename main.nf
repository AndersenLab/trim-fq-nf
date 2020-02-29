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


Channel.fromFilePairs("${params.fastq_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)
  .into { fq_to_QC ; fq_to_trim }




/* 
    ==============================================
    FASTQC on raw data and combine with MultiQC
    ==============================================
*/


process pre_trim_fastqc {

    tag { sampleID }

    publishDir "${params.fastq_path}/${params.fastq_folder}/fastqc", mode: 'copy', pattern: "*.html"

    input:
      set val(sampleID), file(fq1), file(fq2) from fq_to_QC

    output:
       file("*.html") into fastqc_html /* note each item in the channel is a file pair */
       file("*.zip") into fastqc_zip

    """
    fastqc $fq1
    fastqc $fq2
    """
}




process pre_trim_multi_QC {


    publishDir "${params.fastq_path}/${params.fastq_folder}/fastqc", mode: 'copy', pattern: "*.html"

    input:
      file(fastqc_zip) from fastqc_zip.flatten().toSortedList()

    output:
      file("*.html") into pre_multiQC_report
    
    script:

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
      set val(sampleID), file(fq1), file(fq2) from fq_to_trim

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



/* 
    =======================
    combine all trim report
    =======================
*/



process multi_QC {


    publishDir "${params.trimmed_path}/${params.trimmed_folder}/multi_QC", mode: 'copy', pattern: "*.html"

    input:
      file(json) from trim_json.toSortedList()

    output:
      file("*.html") into multiQC_report
    
    script:
      """

      multiqc .

      """
}
