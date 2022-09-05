setwd("/home/amelie/nas/microbiome/res_taxa_scaled/")
rm(list=ls())

res <- fread("res_bayes.txt") %>% 
  dplyr::select(taxa, beta_int, sd_int) %>% 
  mutate(prob=pnorm(0.1, abs(beta_int), sd_int, lower.tail=FALSE), OR=exp(beta_int)) %>% 
  arrange(desc(prob))

res_above05 <- res %>% 
  filter(prob>0.5) #%>% 
  #write.table("res_probabove0.5.txt", row.names = F, quote = F)

res %>% 
  write.table("res_all.txt", row.names = F, quote = F)

res_above05

library(tidyverse)
library(phyloseq)
library(copiome)

load("../dada2_gtdb_2021_11_18_all.RData", verbose = TRUE)

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

taxa_table <- tax_table(phy) %>% as.data.frame()

taxa_above05 <- taxa_table[res_above05$taxa,]
write.table(taxa_above05, "taxa_above05.txt", row.names = T, quote = F)
fread("taxa_above05.txt") %>% 
  write.table("taxa_above05.txt", row.names = F, quote = F)
