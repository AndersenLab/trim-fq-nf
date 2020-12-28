# trim-fq-nf

### Typical use for debugging:

```
nextflow main.nf --debug
```

### Typical use for new fastq:
```
nextflow main.nf --fastq_folder 20180405_fromNUSeq
```

### Parameters
    parameters              description                                   Set/Default
    ==========              ===========                                   ========================
    --debug                 Use --debug to indicate debug mode            (optional)
    --fastq_folder          Name of the raw fastq folder                  (required)
    --raw_path              Path to raw fastq folder                      /projects/b1059/data/fastq/WI/dna/raw
    --processed_path        Path to processed fastq folder (output)       /projects/b1059/data/fastq/WI/dna/processed
    --trim_only             Whether to skip species check and only trim   false
    --genome_sheet          File with fasta locations for species check   genome_sheet.tsv in main folder
    --species_output        Folder name to write species check results    species_check in current folder
    --subsample_read_count  How many reads to use for sepciec check       10000
