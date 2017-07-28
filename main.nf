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

    publishDir params.directory + "/fastqc", mode: 'move'
    
    afterScript "multiqc ${params.directory}/fastqc"
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file(forward), file(reverse) from pre_trim_fastqc
    
    output:
        set file("${forward}_fastqc.zip"), file("${reverse}_fastqc.zip")
        set file("${forward}_fastqc.html"), file("${reverse}_fastqc.html")
    """
        fastqc --noextract --threads 8 ${forward} ${reverse}
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

    publishDir params.out + "/fastqc", mode: 'move'
    
    afterScript "multiqc ${params.out}/fastqc"
    
    cpus 8
    
    tag { dataset_id } 
    
    input:
        set dataset_id, file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") from trim_output
    
    output:
        set dataset_id, file("${dataset_id}_1P_fastqc.zip"), file("${dataset_id}_2P_fastqc.zip"), file("${dataset_id}_1P_fastqc.html"), file("${dataset_id}_2P_fastqc.html") into post_multiqc
    
    """
        fastqc --noextract --threads 8 ${dataset_id}_1P.fq.gz ${dataset_id}_2P.fq.gz
    """
}
