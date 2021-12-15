#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
    - Dan Lu <dan.lu@northwestern.edu>
    - Katie Evans <kathryn.evans@northwestern.edu>
*/

// make sure nextflow version is 20+
if( !nextflow.version.matches('20+') ) {
    println "This workflow requires Nextflow version 20.0 or greater -- You are running version $nextflow.version"
    println "On QUEST, you can use `module load python/anaconda3.6; source activate /projects/b1059/software/conda_envs/nf20_env`"
    exit 1
}

nextflow.preview.dsl=2
// NXF_VER=20.01.0" Require later version of nextflow
//assert System.getenv("NXF_VER") == "20.01.0"

// these aren't used??
include md5sum as md5sum_pre from './md5.module.nf'
include md5sum as md5sum_post from './md5.module.nf'


// read input
if (params.debug) {
    println """

        *** Using debug mode ***

    """
    params.raw_path="${workflow.projectDir}/test_data/raw" 
    params.fastq_folder="20210406_test1"
    params.processed_path="${workflow.projectDir}/test_data/processed"

} else {

    params.raw_path="/projects/b1059/data/transfer/raw"
    params.processed_path="/projects/b1059/data/transfer/processed"
    //params.processed_path="/projects/b1059/data"

}

    params.out="processFQ-${params.fastq_folder}"



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


params.genome_sheet = "${workflow.projectDir}/bin/genome_sheet.tsv"
params.subsample_read_count = "10000"  
md5sum_path = "${params.processed_path}/${params.fastq_folder}/md5sums.txt"
params.R_libpath = "/projects/b1059/software/R_lib_3.6.0"


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

    parameters              description                                   Set/Default
    ==========              ===========                                   ========================
    --debug                 Use --debug to indicate debug mode            ${params.debug}
    --fastq_folder          Name of the raw fastq folder                  ${params.fastq_folder}
    --raw_path              Path to raw fastq folder                      ${params.raw_path}
    --processed_path        Path to processed fastq folder (output)       ${params.processed_path}
    --trim                  Whether to trim fastq                         ${params.trim}
    --check_species         Whether to do species check                   ${params.check_species}
    --genome_sheet          File with fasta locations for species check   ${params.genome_sheet}
    --out                   Folder name to write results                  ${params.out}
    --subsample_read_count  How many reads to use for species check       ${params.subsample_read_count}

    username                                                              ${"whoami".execute().in.text}
    ----------------------------------------------------------------------------------------------
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

"""

out

}


log.info(log_summary())


if (params.help) {
    exit 1
}


println "Running fastp trimming on ${params.raw_path}/${params.fastq_folder}"



workflow { 

    // create sample sheet
    generate_sample_sheet()

    genome_sheet = Channel.fromPath(params.genome_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "genome sheet not found" }
                      .splitCsv(header:true, sep: "\t")

    fq = Channel.fromFilePairs("${params.raw_path}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)
                .concat(Channel.fromFilePairs("${params.raw_path}/${params.fastq_folder}/*_{R1,R2}_*.fastq.gz", flat: true))
                .concat(Channel.fromFilePairs("${params.raw_path}/${params.fastq_folder}/*_{1P,2P}.fq.gz", flat: true))


    // screen species
    if("${params.check_species}" == true) {
    	fq.combine(genome_sheet) | screen_species
        screen_species.out.collect() | multi_QC_species

        // run more species check and generate species-specific sample sheet
        generate_sample_sheet.out
            .combine(multi_QC_species.out) | species_check
    }

    // fastp trim
    if("${params.trim}" == true) {
    	fq | fastp_trim
        fastp_trim.out.fastp_json.collect() | multi_QC_trim
    }


}


/* 
    ==================
    trim raw data
    ==================
*/

process fastp_trim {

    tag { sampleID }

    publishDir "${params.processed_path}/${params.fastq_folder}", mode: 'copy', pattern: "*.fq.gz"

    publishDir "${params.out}/multi_QC", mode: 'copy', pattern: "*_fastp.html"

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
        path("*.stats")

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

    publishDir "${params.out}/multi_QC", mode: 'copy'

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

    publishDir "${params.out}", mode: 'copy'

    input:
        path("*")

    output:
        path("multiqc_data/multiqc_samtools_stats.txt")

    """
        multiqc .
    """
}

/* 
    ================================
    Make sample sheet
    ================================
*/

// this only works for R1_001.fastq.gz right now...

process generate_sample_sheet {
    
    publishDir "${params.out}", mode: 'copy'

    output:
        path("sample_sheet_${params.fastq_folder}_all_temp.tsv")
        

    """
    fq_sheet=`mktemp`
    date=`echo ${params.fastq_folder} | cut -d _ -f 1`
    prefix="${params.raw_path}/${params.fastq_folder}"

    ls \${prefix}/*.gz -1 | xargs -n1 basename | \
    awk -v prefix=\${prefix} -v seq_folder=${params.fastq_folder} -v date=\$date '{
        fq1 = \$1;
        fq2 = \$1;
        gsub("R1_001.fastq.gz", "R2_001.fastq.gz", fq2);
        gsub("R1_001.fastq.gz", "1P.fq.gz", fq1);
        gsub("R2_001.fastq.gz", "2P.fq.gz", fq2);
        split(\$0, a, "_");
        SM = a[1];
        split(\$0, b, "_R");
        ID = b[1];
        gsub("\$", "_", ID);
        gsub("\$", date, ID);
        LB = b[1]
        gsub("\$", "_", LB);
        gsub("\$", date, LB);
        print SM"\\t"ID"\\t"LB"\\t"fq1"\\t"fq2"\\t"seq_folder
    }' | sed -n '1~2p' >> \${fq_sheet}

    if [[ \$(cut -f 2 \${fq_sheet} | sort | uniq -c | grep -v '1 ') ]]; then
        >&2 echo "There are duplicate IDs in the sample sheet. Please review 'inventory.error'"
    fi

    cat \${fq_sheet} | sort | sed '1 i\\strain\\tid\\tlb\\tfq1\\tfq2\\tseq_folder' > sample_sheet_${params.fastq_folder}_all_temp.tsv

    """

}

//         print SM"\\t"ID"\\t"LB"\\t"prefix"/"fq1"\\t"prefix"/"fq2"\\t"seq_folder


/* 
    ================================
    analyze species check
    ================================
*/

process species_check {

    publishDir "${params.out}/sample_sheet/", mode: 'copy', pattern: 'sample_sheet*.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*multiple_libraries.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*master_sheet.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: 'WI_all*.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*_species.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*.html'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*.Rmd'

    input:
        tuple file("sample_sheet"), file("multiqc_samtools_stats")

    output:
        tuple file("*.tsv"), file("*.Rmd"), file("*.html")

    """
        # for some reason, tsv aren't being saved in r markdown, so get around with an R script
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" | cat - ${workflow.projectDir}/bin/species_check.R > species_check.R 
        Rscript --vanilla species_check.R ${params.fastq_folder} ${multiqc_samtools_stats} ${sample_sheet}

        # copy R markdown and insert date and pool
        cat "${workflow.projectDir}/bin/species_check.Rmd" | sed "s/FQ_HOLDER/${params.fastq_folder}/g" > species_check_${params.fastq_folder}.Rmd 

        # add R library path
        echo ".libPaths(c(\\"${params.R_libpath}\\", .libPaths() ))" > .Rprofile

        # make markdown
        Rscript -e "rmarkdown::render('species_check_${params.fastq_folder}.Rmd', knit_root_dir='${workflow.launchDir}')"

    """


}

workflow.onComplete {

    summary = """

    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]

    { Parameters }
    ---------------------------
    --debug                     ${params.debug}
    --fastq_folder              ${params.fastq_folder}
    --raw_path                  ${params.raw_path}
    --processed_path            ${params.processed_path}
    --trim                      ${params.trim}
    --check_species             ${params.check_species}
    --genome_sheet              ${params.genome_sheet}
    --out                       ${params.out}
    --subsample_read_count      ${params.subsample_read_count}

    """

    println summary

    def outlog = new File("${params.out}/log.txt")
    outlog.newWriter().withWriter {
        outlog << param_summary
        outlog << summary
    }


}





