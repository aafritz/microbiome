setwd("/home/amelie/nas/microbiome/")
rm(list=ls())

library(tidyverse)
library(phyloseq)
library(copiome)

load("dada2_gtdb_2021_11_18_all.RData", verbose = TRUE)

phy
# OTU - Operational Taxononic Unit
# New method (DADA2) instead makes ASV - Amplicon Sequence Variant

# First slot = count matrix (otu_table)
taxa_are_rows(phy)
otu_table(phy)[1:5,1:5]

# 2nd slot - sample data
head(sample_data(phy))

# 3rd slot - tax table. Annotated with GTDB https://gtdb.ecogenomic.org/
head(tax_table(phy))

# 4th slot
phy_tree(phy)

# Final slot - the referece sequences
refseq(phy)

# tutorial and references for phyloseq: https://joey711.github.io/phyloseq/preprocess.html
# Overview of functions: https://joey711.github.io/phyloseq/import-data.html

# So which samples do we have?
phy %>% 
  get_variable %>% 
  count(Sampletype, Time)

# subset f1y - 1 year gut.
f1y <- phy %>% subset_samples(Sampletype == "Faeces" & Time == "1y")
f1y <- f1y %>% prune_taxa(taxa_sums(.) > 0, .)


sample_sums(f1y) %>% qplot # = library size
estimate_richness(f1y, measures = "Observed")
qplot(sample_sums(f1y), estimate_richness(f1y, measures = "Observed")$Observed, log = "x")

# relative abundance
f1y_ra <- f1y %>% transform_sample_counts(function(x) x/sum(x))

?make_mradat

# make relative abundance with function of copiome::make_mradat
f1y_mra <- make_mradat(f1y)

ggplot(f1y_mra, aes(prevalence)) + geom_histogram()

f1y_g <- tax_glom(f1y, "Genus")

# make transfer to log + add a number instead of 0
f1y_logrelabu <- f1y_ra %>% otu_df %>% as_tibble %>% 
  mutate_at(vars(everything()), 
            list(function(x) log(x + min(x[x != 0])))) # add lowest nonzero value as pseudocount to prevent -Inf

taxa_over_20pct <- f1y_mra %>% 
  filter(prevalence > .2) %>% 
  pull(tax)

f1y_logrelabu_20pct <- f1y_logrelabu %>% select(all_of(taxa_over_20pct)) %>% 
  bind_cols(abcno = get_variable(f1y, "ABCNO"), .)

f1y_logrelabu_20pct$`6ec6d03fbef9f16e3581ccdc60e7d266` %>% qplot

# "Differential abundance" - benchmarking in your own data https://github.com/Russel88/DAtest/wiki