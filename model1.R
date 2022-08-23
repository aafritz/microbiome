setwd("/home/amelie/nas/microbiome/")
rm(list=ls())

library(tidyverse)
library(phyloseq)
library(copiome)
library(readxl)

load("dada2_gtdb_2021_11_18_all.RData", verbose = TRUE)

# subset f1y - 1 year gut.
f1y <- phy %>% subset_samples(Sampletype == "Faeces" & Time == "1y")
f1y <- f1y %>% prune_taxa(taxa_sums(.) > 0, .)

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
f1y_logrelabu_20pct <- f1y_logrelabu %>% select(all_of(taxa_over_20pct)) %>% 
  bind_cols(abcno = get_variable(f1y, "ABCNO"), .)

## add phenotypes
abcno_pheno <- read_excel("ABC0232_J45_cox_cross.xlsx", sheet = "Data") %>% 
  select(ABCNO, j45_6yr_ever)

## join phenotyptes to taxa table & remove ind with NA
f1y_logrelabu_20pct_pheno <- f1y_logrelabu_20pct %>% 
  left_join(abcno_pheno, by=c("abcno"="ABCNO")) %>% 
  na.omit()

write.table(f1y_logrelabu_20pct_pheno, "f1y_logrelabu_20pct_pheno.txt", row.names = F, quote = F)
