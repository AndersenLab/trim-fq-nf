#!/usr/bin/env nextflow

params.directory = "$PWD/"
params.out = params.directory.replace("raw", "processed")
println params.out
params.threads = 8
println "Running Trimmomatic on " + params.directory
println params.directory + '*_R{1,2}_.*.fastq.gz'

// Fetch fqs; alternative suffixes
Channel.fromFilePairs(params.directory + '*1.fq.gz', flat: true)
        .into { trimmomatic_read_pairs; log_fq }
        
log_fq.subscribe { println it }

// Round up fastq pairs
fastq_suffix.concat( fastq_suffix2, fq_suffix, fq_suffix2 ).into { trimmomatic_read_pairs }

process make_out_dir {
    
    executor 'local'

    """
    mkdir -p ${params.out}
    """
}

process trim {

    publishDir params.out, mode: 'move'

    cpus 8

    input:
        set dataset_id, file(forward), file(reverse) from trimmomatic_read_pairs

    output:
        set file("${dataset_id}_1P.fq.gz"), file("${dataset_id}_2P.fq.gz") into trim_output

    """
    trimmomatic PE -threads ${params.threads} $forward $reverse -baseout ${dataset_id}.fq.gz ILLUMINACLIP:/home/dec211/.linuxbrew/share/trimmomatic/adapters/NexteraPE-PE.fa:2:80:10 MINLEN:45
    rm ${dataset_id}_1U.fq.gz
    rm ${dataset_id}_2U.fq.gz
    """

}

process perform_fq_profile {
    input:
        set f1, f2 from trim_output

    """
    fq profile --fastqc ${params.out}\$(basename ${f1}) ${params.out}/\$(basename ${f2})
    """
}

