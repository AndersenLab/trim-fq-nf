#!/usr/bin/env nextflow

params.directory = "$PWD/"
params.out = params.directory.replace("raw", "processed")
println params.out
params.threads = 8
println "Running Trimmomatic on " + params.directory
println params.directory + '*_{1,2}.fq.gz'

// Fetch fqs; alternative suffixes
Channel.fromFilePairs(params.directory + '*_{1,2}.fq.gz', flat: true)
        .into { pre_trim_fastqc; trimmomatic_read_pairs; log_fq }
        
log_fq.subscribe { println it }

process make_out_dir {
    
    executor 'local'
    
    """
    mkdir -p ${params.out}
    """
}

/*
    Perform fastqc prior to trimming
*/

process pre_trim_fastqc {

    publishDir params.directory + "/fastqc", mode: 'copy'
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file(forward), file(reverse) from pre_trim_fastqc
    
    output:
        set dataset_id, file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_2_fastqc.html"), file("${dataset_id}_2_fastqc.html") into pre_trim_multi_qc
    """
        fastqc --noextract --threads 8 ${forward}
        fastqc --noextract --threads 8 ${reverse}
    """
}

process pre_trim_multi_qc_run {

    publishDir params.directory + "/fastqc", mode: 'copy'

    input:
        set dataset_id, file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_2_fastqc.html"), file("${dataset_id}_2_fastqc.html") from pre_trim_multi_qc.toSortedList()
    output:
        file("multiqc_report.html")
    """
        multiqc .
    """

}

process trim {

    publishDir params.out, mode: 'move'

    cpus 8
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
        set dataset_id, file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") into trim_output

    """
    trimmomatic PE -threads ${params.threads} $forward $reverse -baseout ${dataset_id}.fq.gz ILLUMINACLIP:/home/dec211/.linuxbrew/share/trimmomatic/adapters/NexteraPE-PE.fa:2:80:10 MINLEN:45
    rm ${dataset_id}_1U.fq.gz
    rm ${dataset_id}_2U.fq.gz
    """

}

process post_trim_fastqc {

    publishDir params.out + "/fastqc", mode: 'copy'
    
    validExitStatus 0,2
    
    afterScript "multiqc ${params.out}/fastqc"
    
    cpus 8
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") from trim_output
    
    output:
        set file("${dataset_id}_1P_fastqc.zip"), file("${dataset_id}_2P_fastqc.zip"), file("${dataset_id}_1P_fastqc.html"), file("${dataset_id}_2P_fastqc.html") into post_trim_multi_qc
    
    """
        fastqc --noextract --threads 8 ${dataset_id}_1P.fq.gz
        fastqc --noextract --threads 8 ${dataset_id}_2P.fq.gz
    """
}

process post_trim_multi_qc_run {

    publishDir params.out + "/fastqc", mode: 'copy'

    input:
        set dataset_id, file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_1_fastqc.zip"), file("${dataset_id}_2_fastqc.html"), file("${dataset_id}_2_fastqc.html") from pre_trim_multi_qc.toSortedList()
    
    output:
        file("multiqc_report.html")
        
    """
        multiqc .
    """

}
