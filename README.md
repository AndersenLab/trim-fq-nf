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
   
```

____           .__                     _____                            _____ 
_/  |_ _______ |__|  _____           _/ ____\\  ______           ____  _/ ____\\
\\   __\\\\_  __ \\|  | /     \\   ______ \\   __\\  / ____/  ______  /    \\ \\   __\\ 
 |  |   |  | \\/|  ||  Y Y  \\ /_____/  |  |   < <_|  | /_____/ |   |  \\ |  |   
 |__|   |__|   |__||__|_|  /          |__|    \\__   |         |___|  / |__|   
                         \\/                      |__|              \\/         
                                                                              

    parameters              description                                   Set/Default
    ==========              ===========                                   ========================
    --debug                 Use --debug to indicate debug mode            ${params.debug}
    --fastq_folder          Name of the raw fastq folder                  ${params.fastq_folder}
    --raw_path              Path to raw fastq folder                      ${params.raw_path}
    --processed_path        Path to processed fastq folder (output)       ${params.processed_path}
    --genome_sheet          File with fasta locations for species check   ${params.genome_sheet}
    --out                   Folder name to write results                  ${params.out}
    --subsample_read_count  How many reads to use for species check       ${params.subsample_read_count}
    
```


## Software requirements

* Nextflow v20.01+ (see the dry guide on Nextflow [here](quest-nextflow.md) or the Nextflow documentation [here](https://www.nextflow.io/docs/latest/getstarted.html)). On QUEST, you can access this version by loading the `nf20` conda environment prior to running the pipeline command:

```
module load python/anaconda3.6
source activate /projects/b1059/software/conda_envs/nf20_env
```

* Docker - this pipeline uses three separate docker images: 1) `andersenlab/trim-fq` and 2) `andersenlab/multiqc` which are both generated as part of this repo (`env/trim-fq.Dockerfile` for example). Check out the [dry guide](http://andersenlab.org/dry-guide/2021-12-01/pipeline-docker/) for more info. This pipeline also uses the `andersenlab/r_packages` container for all R work, which is hosted and generated [separately](https://github.com/AndersenLab/dockerfile/tree/master/r_packages).
    - If you are on QUEST, you can load docker (actually singularity) with:

```
module load singularity
```

**Note: As of 2022-01-01, the conda environments that used to host this pipeline will no longer be maintained**

# Usage

## Testing the pipeline on QUEST

*This command uses a test dataset*

```
nextflow run main.nf --debug
```

## Running the pipeline on QUEST

```
nextflow run main.nf --fastq_folder <name_of_folder>
```

# Profiles

## -profile standard (Default)

If no profile is designated, the default profile will run both fastq trimming AND species check

## -profile trim_only

Use this profile to only trim fastq files and not perform species check.

## -profile sp_check_only

Use this profile to only run species check and not fastq trimming. This is useful for running species checks on previously trimmed fastqs.

# Parameters

## --debug

You should use `--debug true` for testing/debugging purposes. This will run the debug test set (located in the `test_data/raw` folder).

For example:

```
nextflow run main.nf --debug -resume
```

Using `--debug` will automatically set the fastq_folder to `test_data/raw/20210406_test1`

## --fastq_folder

This should be the name of the folder containing all fastq files located at `/projects/b1059/data/transfer/raw/`. As long as there are no overlapping file names (be sure to check this first), you can combine multiple pools sequenced at the same time into one larger folder at this step.

### --raw_path (optional)

The path to the `fastq_folder` if not default (`/projects/b1059/data/transfer/raw/`)

### --processed_path (optional)

The path to output folder if not default (`/projects/b1059/data/transfer/processed/`)


### --genome_sheet (optional)

Path to a tsv file listing project IDs for species. Default is located in `bin/genome_sheet.tsv`

### --out (optional)

Name of output folder with results. Default is "processFQ-{fastq_folder}"

### --subsample_read_count (optional)

How many reads to use for species check. Default = 10,000


# Output

```
├── b1059/data/transfer/processed/
│   └── {strain}_{library}_{read}.fq.gz
- - - - - - - - - - - - - - - - - - - - - - - - - - - -
├── multi_QC
│   ├── multiqc_data
│   │   ├── *.txt
│   │   └── multiqc_data.json
│   ├── multiqc_report.html
│   └── {strain}_{library}_{read}_fastp.html
├── multiqc_data
│   └── multiqc_samtools_stats.txt
├── sample_sheet
│   ├── sample_sheet_{species}_{date}_ALL.tsv
│   └── sample_sheet_{species}_{date}_NEW.tsv
├── species_check
│   ├── species_check_{fastq_folder}.html
│   ├── {library}_multiple_librarires.tsv
│   ├── {library}_strains_most_likely_species.tsv
│   ├── {library}_strains_not_in_master_sheet.tsv
│   ├── {library}_strains_possibly_diff_species.tsv
│   └── WI_all_{date}.tsv
├── sample_sheet_{fastq_folder}_all_temp.tsv
└── log.txt
```

The resulting trimmed FASTQs will be output in the `b1059/data/transfer/processed` directory. The rest of the output files and reports will be generated in a new folder in the directory in which you ran the nextflow pipeline, labeled by `processFQ-{fastq_folder}`.

__MultiQC__

* `multi_QC/{strain}_{library}_fastp.html` - fastp html report detailing trimming and quality

* `multi_QC/multiqc_report.html` - aggregate multi QC report for all strains pre and post trimming

* `multi_QC/multiqc_data/*.txt or .json` - files used to make previous reports


__Species check__

If a species check is run, the `multiqc_data/multiqc_samtools_stats.txt` will contain the results of reads mapped to each species. Furthermore, several reports and sample sheets will be generated:
 
* `species_check/species_check_{date}_{library}.html` is an HTML report showing how many strains have issues (not in master sheet, possibly different species, etc.)

* `species_check/{library}_multiple_libraries.tsv` - strains sequenced in multiple libraries

* `species_check/{library}_strains_most_likely_species.tsv` - list of all strains in library labeled by (1) the species in record and (2) most likely species by sequencing

* `species_check/{library}_strains_not_in_master.tsv` - list of strains not found in Robyn's wild isolate master sheets for CE, CB, or CT

* `species_check/{library}_strains_possibly_diff_species.tsv` - list of strains whose species in record does not match the likely species by sequencing

* `species_check/WI_all_{date}.tsv` - copy of all strains (CE, CB, and CT) and species designation in record

* `sample_sheet_{date}_{library}_all_temp.tsv` - temporary sample sheet with all strains for all species combined. DO NOT USE THIS FOR ALIGNMENT.

__Sample sheets__

If a species check is run, the `species_check/sample_sheet` folder will also contain 6 sample sheets to be used for alignment:

* `sample_sheet/sample_sheet_{species}_{date}_ALL.tsv` - sample sheet for `alignment-nf` using ALL strains of a particular species (i.e. c_elegans). This is useful for species we have not performed any alignments for or when we update the reference genome and need to re-align all strains.

* `sample_sheet/sample_sheet_{species}_{date}_NEW.tsv` - sample sheet for `alignment-nf` using all fastq from any library for ONLY strains sequenced in this particular library of a particular species (i.e. c_elegans, RET63). This is useful when the reference genome does not change and there is no need to re-align thousands of strains to save on computational power.