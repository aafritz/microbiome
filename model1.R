setwd("/home/amelie/nas/microbiome/")
rm(list=ls())

library(tidyverse)
library(phyloseq)
library(copiome)
library(readxl)
library(zCompositions)
library(compositions)

load("dada2_gtdb_2021_11_18_all.RData", verbose = TRUE)

# subset f1y - 1 year gut.
f1y_0 <- phy %>% subset_samples(Sampletype == "Faeces" & Time == "1y")
f1y <- f1y_0 %>% prune_taxa(taxa_sums(.) > 0, .)

# taxa_sums(): returns the total number of individuals observed from each species/taxa/OTU
# sums up the count of all ind of each taxa
# prune_taxa(taxa that should be removed, data set from which they should be removed)

# relative abundance
f1y_ra <- f1y %>% transform_sample_counts(function(x) x/sum(x))

# make relative abundance with function of copiome::make_mradat
f1y_mra <- make_mradat(f1y)

# make transfer to log + add a number instead of 0
f1y_logrelabu <- f1y_ra %>% otu_df %>% as_tibble %>% 
  mutate_at(vars(everything()), 
            list(function(x) log(x + min(x[x != 0])))) # add lowest nonzero value as pseudocount to prevent -Inf

# remove all taxa below 0.2
taxa_over_20pct <- f1y_mra %>% 
  filter(prevalence > .2) %>% 
  pull(tax)

# add ABC numbers
f1y_logrelabu_20pct <- f1y_logrelabu %>% dplyr::select(all_of(taxa_over_20pct)) %>% 
  bind_cols(abcno = get_variable(f1y, "ABCNO"), .)

## add phenotypes
abcno_pheno <- read_excel("ABC0232_J45_cox_cross.xlsx", sheet = "Data") %>% 
  dplyr::select(ABCNO, j45_6yr_ever)

## join phenotyptes to taxa table & remove ind with NA
f1y_logrelabu_20pct_pheno <- abcno_pheno %>% 
  left_join(f1y_logrelabu_20pct, by=c("ABCNO"="abcno")) %>% 
  na.omit()

colnames(f1y_logrelabu_20pct_pheno)[2] <- "pheno"

## add gender
sex <- read_excel("ABC0202_Birth_Sex_Race.xlsx", sheet="Data") %>% 
  dplyr::select(ABCNO, SEX) %>% 
  mutate(SEX=ifelse(SEX=="Male", 0, 1))

f1y_logrelabu_20pct_pheno_sex <- sex %>% 
  left_join(f1y_logrelabu_20pct_pheno, by=c("ABCNO"="ABCNO"))


write.table(f1y_logrelabu_20pct_pheno_sex, "f1y_logrelabu_20pct_pheno_sex.txt", row.names = F, quote = F)
