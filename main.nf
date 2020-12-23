#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
    - Dan Lu <dan.lu@northwestern.edu>
*/

nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"




include md5sum as md5sum_pre from './md5.module.nf'
include md5sum as md5sum_post from './md5.module.nf'


// read input
if (params.debug) {
    println """

        *** Using debug mode ***

    """
    params.raw_path="/projects/b1059/workflows/trim-fq-nf/test_data/raw" 
    params.fastq_folder="fq_run"
    params.processed_path="/projects/b1059/workflows/trim-fq-nf/test_data/processed"

} else {

    params.raw_path="/projects/b1059/data/fastq/WI/dna/raw"
    params.processed_path="/projects/b1059/data/fastq/WI/dna/processed"


}


// required inputs
if (params.fastq_folder == null) {
    if (params.help) {
    } else {
        println """

        Please specify fastq folder name with --fastq_folder

        """
        exit 1
    }
}


params.genome_sheet = "${workflow.projectDir}/genome_sheet.tsv"
params.species_output = "species_check"   // default is to write species output to current folder
params.subsample_read_count = "10000"  // 
md5sum_path = "${params.processed_path}/${params.fastq_folder}/md5sums.txt"



def log_summary() {

  out='''
                                                                              
____           .__                     _____                            _____ 
_/  |_ _______ |__|  _____           _/ ____\\  ______           ____  _/ ____\\
\\   __\\\\_  __ \\|  | /     \\   ______ \\   __\\  / ____/  ______  /    \\ \\   __\\ 
 |  |   |  | \\/|  ||  Y Y  \\ /_____/  |  |   < <_|  | /_____/ |   |  \\ |  |   
 |__|   |__|   |__||__|_|  /          |__|    \\__   |         |___|  / |__|   
                         \\/                      |__|              \\/         
                                                                              

''' + """
To run the pipeline:

nextflow main.nf --debug
nextflow main.nf --fastq_folder 20180405_fromNUSeq

    parameters              description                                  Set/Default
    ==========              ===========                                  ========================
    --debug                 Use --debug to indicate debug mode           ${params.debug}
    --fastq_folder          Name of the raw fastq folder                 ${params.fastq_folder}
    --raw_path              Path to raw fastq folder                     ${params.raw_path}
    --processed_path        Path to processed fastq folder (output)      ${params.processed_path}
    --genome_sheet          File with fasta locations for species check  ${params.genome_sheet}
    --species_output        Folder name to write species check results   ${params.species_output}
    --subsample_read_count  How many reads to use for sepciec check      ${params.subsample_read_count}

    username                                                             ${"whoami".execute().in.text}

"""

out

}


log.info(log_summary())


if (params.help) {
    exit 1
}


println "Running fastp trimming on ${params.raw_path}/${params.fastq_folder}"



workflow { 

    genome_sheet = Channel.fromPath(params.genome_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "genome sheet not found" }
                      .splitCsv(header:true, sep: "\t")

    fq = Channel.fromFilePairs("${params.raw_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)

    fq | fastp_trim
    fastp_trim.out.fastp_json.collect() | multi_QC_trim

    fq.combine(genome_sheet) | screen_species
    screen_species.out.collect() | multi_QC_species



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
    ==================
    screen species
    ==================
*/



process screen_species {

    conda "/projects/b1059/software/conda_envs/alignment-nf_env"

    input:
        tuple val(sampleID), path(fq1), path(fq2), genome_row

    output:
        tuple path("*.stats")

    """
        zcat ${fq1} | head -n ${params.subsample_read_count} > subset_R1.fq
        zcat ${fq2} | head -n ${params.subsample_read_count} > subset_R2.fq
        bwa mem -t ${task.cpus} ${genome_row.genome} subset_R1.fq subset_R2.fq > out.sam

        samtools sort -o sorted.bam out.sam
        samtools index sorted.bam

        picard MarkDuplicates I=sorted.bam \\
                              O=rm_dup.bam \\
                              M=${sampleID}_xxx_${genome_row.species}.duplicates.txt \\
                              VALIDATION_STRINGENCY=LENIENT \\
                              REMOVE_DUPLICATES=true

        samtools stats rm_dup.bam > ${sampleID}_xxx_${genome_row.species}.stats

        rm *.fq
        rm *.bam
        rm *.sam
    """
}



/* 
    =======================
    combine all trim report
    =======================
*/

process multi_QC_trim {

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


/* 
    ================================
    combine all species check report
    ================================
*/


process multi_QC_species {

    publishDir "${params.species_output}", mode: 'copy'

    input:
        path("*")

    output:
        path("multiqc_data/multiqc_samtools_stats.txt")

    """
        multiqc .
    """
}
