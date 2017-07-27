# trimmomatic-nf

The trimmomatic workflow should be used to initially process sequence data. In general, it should be used on
wild isolate sequence data, however - it should not be used on low-coverage sequencing (as used to determine NIL / RIL genotypes). 

`/projects/b1059/data/fastq/WI/dna/raw/<folder_name>`

### Usage

```
# cd to directory of fastqs
nextflow run Andersenlab/trimmomatic-nf
```

### .fq.gz

Fastqs should end with `.fq.gz`; To rename fastqs that end in `fastq.gz` use:

```
rename --dry-run --subst .fastq.gz .fq.gz *.fastq.gz
```

Verify that the results are correct and remove `dry-run` when appropriate.
