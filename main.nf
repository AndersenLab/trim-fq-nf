#!/usr/bin/env nextflow 
/*
    Andersen Lab C. elegans Trimming Pipeline
    Authors:
    - Daniel Cook <danielecook@gmail.com>
    - Dan Lu <dan.lu@northwestern.edu>
    - Katie Evans <kathryn.evans@northwestern.edu>
    - Mike Sauria <mike.sauria@jhu.edu>
*/

// make sure nextflow version is 23+
if( !nextflow.version.matches('23+') ) {
    println "This workflow requires Nextflow version 23.0 or greater -- You are running version $nextflow.version"
    println "On Rockfish, you can use `module load python/anaconda3.6; source activate /vast/eande106/software/conda_envs/nf23_env`"
    exit 1
}

nextflow.enable.dsl=2

// these aren't used??
//include { md5sum as md5sum_pre } from './md5.module.nf'
//include { md5sum as md5sum_post } from './md5.module.nf'


// read input
if (params.debug) {
    println """

        *** Using debug mode ***

    """
    params.raw_path_final="${workflow.projectDir}/test_data/raw" 
    params.fastq_folder="MMDDYYYY_testrun"
    params.processed_path_final="${workflow.launchDir}/processed"

} else {
    params.raw_path_final = params.raw_path
    params.processed_path_final = params.processed_path
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

if ((params.raw_path_final == null) | (params.processed_path_final == null) | (params.genome_path == null)) {
    if (params.help) {
    } else {
        println """

        If running locally, please specify raw, processed, and genome paths manually with --raw_path, --processed_path, and --genome_path

        """
        exit 1
    }
}


params.genome_sheet = "${workflow.projectDir}/bin/genome_sheet.tsv"
params.subsample_read_count = "10000"  
//md5sum_path = "${params.processed_path_final}/${params.fastq_folder_final}/md5sums.txt"


def log_summary() {

  out='''
                                                                              
____           .__                     _____                            _____ 
_/  |_ _______ |__|  _____           _/ ____\\  ______           ____  _/ ____\\
\\   __\\\\_  __ \\|  | /     \\   ______ \\   __\\  / ____/  ______  /    \\ \\   __\\ 
 |  |   |  | \\/|  ||  Y Y  \\ /_____/  |  |   | |_|  | /_____/ |   |  \\ |  |   
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
    --raw_path              Path to raw fastq folder                      ${params.raw_path_final}
    --processed_path        Path to processed fastq folder (output)       ${params.processed_path_final}
    --trim                  Whether to trim fastq                         ${params.trim}
    --check_species         Whether to do species check                   ${params.check_species}
    --genome_sheet          File with fasta locations for species check   ${params.genome_sheet}
    --genome_path           Path to data in genome_sheet                  ${params.genome_path}
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


println "Running fastp trimming on ${params.raw_path_final}/${params.fastq_folder}"



workflow { 

    // create sample sheet
    Channel.fromPath("${params.raw_path_final}/${params.fastq_folder}") | generate_sample_sheet
    
    genome_sheet = Channel.fromPath(params.genome_sheet, checkIfExists: true)
                      .ifEmpty { exit 1, "genome sheet not found" }
                      .splitCsv(header:true, sep: "\t")
    
    genome_names = genome_sheet
        .map{ it: it.species }
        .collect()

    genome_paths = genome_sheet
        .map{ it: it.genome }
        .collect()

    fq = Channel.fromFilePairs("${params.raw_path_final}/${params.fastq_folder}/*_{1,2}.fq.gz", flat: true)
                .concat(Channel.fromFilePairs("${params.raw_path_final}/${params.fastq_folder}/*_{R1,R2}_*.fastq.gz", flat: true))
                .concat(Channel.fromFilePairs("${params.raw_path_final}/${params.fastq_folder}/*_{1P,2P}.fq.gz", flat: true))

    // screen species
    if("${params.check_species}" == true) {
    	screen_species( fq,
                        genome_names,
                        genome_paths,
                        Channel.fromPath(params.genome_path).first() )
        screen_species.out.collect() | multi_QC_species

        // run more species check and generate species-specific sample sheet
        generate_sample_sheet.out
            .combine(multi_QC_species.out)
            .combine(Channel.fromPath("${workflow.projectDir}/bin/species_check.Rmd")) | species_check
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

    publishDir "${params.processed_path_final}/${params.fastq_folder}", mode: 'copy', pattern: "*.fq.gz"

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
    array 100

    input:
        tuple val(sampleID), path(fq1), path(fq2)
        val genome_names
        val genome_indices
        path genome_dir

    output:
        path("*.stats")

    """
    NAMES=(`echo "${genome_names}" | sed "s/,//g" | sed "s/\\[//" | sed "s/\\]//"`)
    INDICES=(`echo "${genome_indices}" | sed "s/,//g" | sed "s/\\[//" | sed "s/\\]//"`)
    N=\${#NAMES[*]}

    zcat ${fq1} | head -n ${params.subsample_read_count} > subset_R1.fq
    zcat ${fq2} | head -n ${params.subsample_read_count} > subset_R2.fq
    
    for I in \$(seq 0 1 \$(expr \${N} - 1)); do
        SPECIES=\${NAMES[\${I}]}
        INDEX=${genome_dir}/\${INDICES[\${I}]}
        bwa mem -t ${task.cpus} \${INDEX} subset_R1.fq subset_R2.fq > out.sam

        samtools sort -o sorted.bam out.sam
        samtools index sorted.bam

        picard MarkDuplicates I=sorted.bam \\
                                O=rm_dup.bam \\
                                M=${sampleID}_xxx_\${SPECIES}.duplicates.txt \\
                                VALIDATION_STRINGENCY=LENIENT \\
                                REMOVE_DUPLICATES=true

        samtools stats rm_dup.bam > ${sampleID}_xxx_\${SPECIES}.stats
    done

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

    // this process uses a different container than the others
    container 'andersenlab/multiqc:2022030115492310c8da'

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

    // this process uses a different container than the others
    container 'andersenlab/multiqc:2022030115492310c8da'

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

    input:
        path(fq_folder)

    output:
        path("sample_sheet_${params.fastq_folder}_all_temp.tsv")
        

    """
    fq_sheet=`mktemp`
    date=`echo ${params.fastq_folder} | cut -d _ -f 1`
    prefix="${params.raw_path_final}/${params.fastq_folder}"

    ls ${fq_folder}/*.gz -1 | xargs -n1 basename | \\
    awk -v prefix=\${prefix} -v seq_folder=${params.fastq_folder} -v date=\${date} \\
    'BEGIN{OFS="\\t"}
    {
        if ((\$1 ~ /.*R1_[0-9]+\\.fastq\\.gz\$/)){
            fq1 = \$1;
            gsub(/R1_[0-9]+\\.fastq/,"1P.fq", fq1);
            fq2 = fq1;
            gsub(/1P.fq/, "2P.fq", fq2);
            valid = "True";
        } else if ((\$1 ~ /.*_1\\.fq\\.gz\$/)){
            fq1 = \$1;
            gsub(/1\\.fq/,"1P.fq", fq1);
            fq2 = fq1;
            gsub(/1P\\.fq/, "2P.fq", fq2);
            valid = "True";
        } else if ((\$1 ~ /.*1P\\.fq\\.gz\$/)){
            fq1 = \$1;
            fq2 = fq1;
            gsub(/1P\\.fq/, "2P.fq", fq2);
            valid = "True";
        } else
            valid = "False";
        if (valid == "True"){
            split(fq1, a, "_");
            SM = a[1];
            split(fq1, b, "_1P");
            ID = b[1];
            gsub("\$", "_", ID);
            gsub("\$", date, ID);
            LB = b[1];
            gsub("\$", "_", LB);
            gsub("\$", date, LB);
            print SM,ID,LB,fq1,fq2;
        }
    }' >> \${fq_sheet}
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

    // use r_packages container
    container 'andersenlab/r_packages:v0.7'

    publishDir "${params.out}/sample_sheet/", mode: 'copy', pattern: 'sample_sheet*.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*multiple_libraries.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*master_sheet.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: 'WI_all*.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*_species.tsv'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*.html'
    publishDir "${params.out}/species_check/", mode: 'copy', pattern: '*.Rmd'

    input:
        tuple file(sample_sheet), file(multiqc_samtools_stats), file(report)

    output:
        tuple file("*.tsv"), file("*.Rmd"), file("*.html")

    """
        # for some reason, tsv aren't being saved in r markdown, so get around with an R script
        Rscript --vanilla ${workflow.projectDir}/bin/species_check.R ${params.fastq_folder} ${multiqc_samtools_stats} ${sample_sheet}

        # copy R markdown and insert date and pool
        cat "${report}" | sed "s/FQ_HOLDER/${params.fastq_folder}/g" > species_check_${params.fastq_folder}.Rmd 

        # make markdown
        Rscript -e "rmarkdown::render('species_check_${params.fastq_folder}.Rmd')"

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
    --raw_path                  ${params.raw_path_final}
    --processed_path            ${params.processed_path_final}
    --trim                      ${params.trim}
    --check_species             ${params.check_species}
    --genome_sheet              ${params.genome_sheet}
    --genome_path               ${params.genome_path}
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





