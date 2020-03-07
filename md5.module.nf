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
      md5_comm=md5sum
    else
      md5_comm=md5
    fi;
    \${md5_comm} ${fq1} > md5.txt
    \${md5_comm} ${fq2} >> md5.txt
    """
}
