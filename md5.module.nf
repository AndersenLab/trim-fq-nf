// A process for generating an md5sum hash

process md5sum {
    // Generate an md5sum of each FASTQ

    input:
      tuple val(sampleID), path(fq1), path(fq2) 
    
    output:
      path "md5.txt"
    
    """
    md5sum ${fq1} > md5.txt
    md5sum ${fq2} >> md5.txt
    """
}