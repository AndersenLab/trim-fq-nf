#!/usr/bin/env Rscript

library(tidyverse)
library(fuzzyjoin)
library(gsheet)

args = commandArgs(trailingOnly=TRUE)

#######################
#  ARGUMENTS
#######################

fq_folder <- args[1]
date_label <- stringr::str_split_fixed(fq_folder, "_", 2)[,1]
RET <- stringr::str_split_fixed(fq_folder, "_", 2)[,2]

#######################
#  LOAD WI DATA
#######################
# Load sample sheet and append species. Note the species name should be the same as what the species_check used
ce_WI <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/10x-CcKNCl80F9hMcrGWC4fhP_cbekSzi5_IYBY2UqCc/edit?usp=sharing") %>%
    dplyr::select(strain) %>%
    dplyr::mutate(species_in_record="c_elegans")

cb_WI <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1IJHMLwuaxS_sEO31TyK5NLxPX7_qSd0bHNKverAv8-0/edit?usp=sharing") %>%
    dplyr::select(strain) %>%
    dplyr::mutate(species_in_record="c_briggsae")

ct_WI <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1mqXOlUX7UeiPBe8jfAwFZnqlzhb7X-eKGK_TydT7Gx4/edit?usp=sharing") %>%
    dplyr::select(strain) %>%
    dplyr::mutate(species_in_record="c_tropicalis")

wi_all = bind_rows(ce_WI, cb_WI, ct_WI)

# load Robyn's sequencing sheet
seq_sheet <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/edit#gid=484600292") %>%
    dplyr::select(sample, species) %>% 
    tidyr::separate(sample, c("strain","pool_in_sheet"), sep = "_")

# keep a record
write.table(wi_all, glue::glue("WI_all_{date_label}.tsv"), quote=F, col.names = T, row.names = F, sep="\t")

#######################
#  LOAD FQ STATS
#######################

# read in multiqc data
sc_max <- read.delim(args[2], stringsAsFactors=FALSE) %>% 
    dplyr::distinct() %>%
    dplyr::filter(Sample != "Sample") %>%
    tidyr::separate(Sample, c("sample","likely_species"), sep="_xxx_") %>% 
    dplyr::mutate(strain = stringr::str_split_fixed(sample, "_", 2)[,1]) %>%
    dplyr::mutate(bases_mapped_percent = as.numeric(bases_mapped_.cigar.)/as.numeric(total_length) * 100) %>% 
    dplyr::select(strain, sample, likely_species, bases_mapped_percent) %>%
    dplyr::group_by(sample) %>% 
    dplyr::top_n(n = 1, wt = bases_mapped_percent) %>% 
    dplyr::left_join(wi_all) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(species_in_record = ifelse(is.na(species_in_record), "not_in_record", species_in_record))

# split out those that are not in master sheet
sc_max_in <- sc_max %>%
    dplyr::filter(species_in_record != "not_in_record")

# Use Robyn's sheet for the species, but need to re-format species names
sc_max_not_in <- sc_max %>%
    dplyr::filter(species_in_record == "not_in_record") %>%
    dplyr::left_join(seq_sheet, by="strain") %>%
    dplyr::select(-species_in_record, -pool_in_sheet) %>%
    dplyr::rename(species_in_record=species) %>% unique() %>%
    dplyr::mutate(species_in_record = str_replace(species_in_record, "C. nigoni", "c_nigoni")) %>%
    dplyr::mutate(species_in_record = str_replace(species_in_record, "C. briggsae", "c_briggsae")) %>% 
    dplyr::mutate(species_in_record = str_replace(species_in_record, "C. elegans", "c_elegans")) %>% 
    dplyr::mutate(species_in_record = str_replace(species_in_record, "C. tropicalis", "c_tropicalis"))

# most likely species
sc_max <- dplyr::bind_rows(sc_max_not_in, sc_max_in)
write.table(sc_max, paste0(RET, "_strains_most_likely_species.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

#######################
#  PROBLEM STRAINS
#######################

# mislabeled strains
problem_strains <- sc_max %>%
    dplyr::filter(likely_species != species_in_record)

# write out problematic strains to put on Robyn's Mislabeled species sheet
problem_strains %>%
    dplyr::filter(species_in_record =="not_in_record") %>%
    write.table( paste0(RET, "_strains_not_in_master_sheet.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

problem_strains %>%
    dplyr::filter(species_in_record !="not_in_record") %>% 
    write.table( paste0(RET, "_strains_possibly_diff_species.tsv"), col.names=T, row.names=F, sep="\t", quote=F)

# write out strains with more than 1 libraries
sc_max %>%
    dplyr::filter(duplicated(strain) | duplicated(strain, fromLast = T)) %>% 
    dplyr::arrange(strain) %>% 
    write.table(paste0(RET, "_multiple_libraries.tsv"), sep="\t", quote=F, col.names = T, row.names=F)

#######################
#  CREATE SAMPLE SHEET
#######################

# add max-mapped species to sample sheet
ss <- read.delim(args[3], stringsAsFactors=FALSE)

sc_max <- sc_max %>%
    dplyr::select(sample, likely_species, bases_mapped_percent, species_in_record) 

all <- fuzzyjoin::fuzzy_inner_join(ss, sc_max, by = c("fq1" = "sample"), match_fun = str_detect) 

# C. elegans
c_elegans_new <- all %>%
    dplyr::filter(likely_species=="c_elegans" & species_in_record=="c_elegans") %>% 
    dplyr::select(-likely_species, -bases_mapped_percent, -species_in_record, -sample) %>% 
    dplyr::arrange(seq_folder, strain)
    # write.table(glue::glue( "fq_sheet_for_seq_sheet_CE_{date_label}.tsv"), sep="\t", quote=F, col.names = T, row.names = F)

# C. tropicalis
c_tropicalis_new = all %>%
    dplyr::filter(likely_species=="c_tropicalis" & species_in_record=="c_tropicalis") %>% 
    dplyr::select(-likely_species, -bases_mapped_percent, -species_in_record, -sample) %>% 
    dplyr::arrange(seq_folder, strain) 
    # write.table(glue::glue( "fq_sheet_for_seq_sheet_CT_{date_label}.tsv"), sep="\t", quote=F, col.names = T, row.names = F)

# C. briggsae
c_briggsae_new = all %>%
    dplyr::filter(likely_species=="c_briggsae" & species_in_record=="c_briggsae") %>% 
    dplyr::select(-likely_species, -bases_mapped_percent, -species_in_record, -sample) %>% 
    dplyr::arrange(seq_folder, strain)
    # write.table(glue::glue( "fq_sheet_for_seq_sheet_CB_{date_label}.tsv"), sep="\t", quote=F, col.names = T, row.names = F)

#######################
#  CREATE FULL SAMPLE SHEET
#######################
# load sheet containing all fq before this pool
fq_elegans <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/edit#gid=719017664") %>%
    dplyr::select(strain:seq_folder) %>%
    dplyr::mutate(species = "c_elegans")
fq_briggsae <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/edit#gid=872003994") %>%
    dplyr::select(strain:seq_folder) %>%
    dplyr::mutate(species = "c_briggsae")
fq_tropicalis <- gsheet::gsheet2tbl("https://docs.google.com/spreadsheets/d/1CpSpzU1p-WtGKIMBK99DL5AeZb-A8QrHPuLkM_fAuEY/edit#gid=1811021622") %>%
    dplyr::select(strain:seq_folder) %>%
    dplyr::mutate(species = "c_tropicalis")

fq_all <- fq_elegans %>%
    dplyr::bind_rows(fq_briggsae, fq_tropicalis)

# keep a record
write.table(fq_all, glue::glue("{date_label}_fq_sheet_all.tsv"), sep="\t", quote=F, col.names = T, row.names = F)

# write out sample sheet for alignment
for(ss in unique(fq_all$species)) {
    # get all fastq for strains in this pool
    fq_new <- get(glue::glue("{ss}_new"))
    fq_new1 <- fq_all %>%
        dplyr::filter(species == ss) %>%
        dplyr::select(-species) %>%
        dplyr::semi_join(fq_new, by="strain")
    fq_new2 <- fq_new %>%
        dplyr::bind_rows(fq_new1) %>% 
        unique() %>% 
        dplyr::arrange(strain) # unique here in case the new fq are already appended
    
    # this includes fq from previous pools of the same strains in this pool
    write.table(fq_new2, glue::glue("sample_sheet_{ss}_{date_label}_NEW.tsv"), sep="\t", quote=F, col.names = T, row.names = F)

    # this is all fq for species
    total <- fq_all %>%
        dplyr::filter(species == ss) %>%
        dplyr::select(-species) %>%
        dplyr::bind_rows(fq_new)
    
    write.table(total, glue::glue("sample_sheet_{ss}_{date_label}_ALL.tsv"), sep="\t", quote=F, col.names = T, row.names = F)


    #######################
    #  APPEND TO GOOGLESHEETS
    #######################   
    # NOT HERE YET
}
