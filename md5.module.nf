// A process for generating an md5sum hash

process md5sum {
    // Generate an md5sum of each FASTQ

    input:
      tuple val(sampleID), path(fq1), path(fq2) 
    
    output:
      path "md5.txt"
    
    """
    # Allow this command to work on unix/macos
    if command -v md5sum; then
      alias md5=md5sum
    fi;
    md5 ${fq1} > md5.txt
    md5 ${fq2} >> md5.txt
    """
}