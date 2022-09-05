setwd("/home/amelie/nas/microbiome/")
rm(list=ls())

library(tidyverse)
library(phyloseq)
library(copiome)
library(readxl)
library(zCompositions)
library(compositions)
library(qqman)

load("dada2_gtdb_2021_11_18_all.RData", verbose = TRUE)

# subset f1y - 1 year gut.
f1y_0 <- phy %>% subset_samples(Sampletype == "Faeces" & Time == "1y")
f1y <- f1y_0 %>% prune_taxa(taxa_sums(.) > 0, .)

# taxa_sums(): returns the total number of individuals observed from each species/taxa/OTU
# sums up the count of all ind of each taxa
# prune_taxa(taxa that should be removed, data set from which they should be removed)

#### AT GENUS LEVEL ####
f1y_g <- tax_glom(f1y, "Genus") %>% 
  otu_df()

f1y_g.GBM <- cmultRepl(otu_table(f1y_g))
write.table(f1y_g.GBM, "f1y_g.GBM.txt", quote = F)
f1y_g.GBM <- fread("f1y_g.GBM.txt") %>% 
  dplyr::rename(taxa=V1)

# centered log transform
f1y_g.GBM.clr <- clr(f1y_g.GBM[,-1])
f1y_g.GBM.clr.df <- as.data.frame(f1y_g.GBM.clr)
f1y_g.GBM.clr2 <- cbind(f1y_g.GBM[,1], f1y_g.GBM.clr.df)

# select all genera > .2
f1y_g_mra <- make_mradat(f1y_g)
taxa_over_20pct_g <- f1y_g_mra %>% 
  filter(prevalence > .2) %>% 
  pull(tax)

f1y_g.GBM.clr2_20pct_0 <- subset(f1y_g.GBM.clr2, !(taxa %in% taxa_over_20pct_g))
  
f1y_g.GBM.clr2_20pct <- f1y_g.GBM.clr2_20pct_0[,-1] %>% 
  bind_cols(abcno = get_variable(f1y, "ABCNO"), .)

f1y_g.GBM.clr2 %>% subset(taxa %in% taxa_over_20pct_g)
#### AT SPECIES LEVEL ####
# relative abundance with zcomposition package
f1y.GBM <- cmultRepl(otu_table(f1y))

# centered log ratio transformation
f1y.GBM.clr <- clr(f1y.GBM[,-1])
f1y.GBM.clr.df <- as.data.frame(f1y.GBM.clr)
f1y.GBM.clr2 <- cbind(f1y.GBM[,1], f1y.GBM.clr.df) %>% 
  dplyr::rename(taxa=V1)

write.table(f1y.GBM.clr2, "f1y.GBM.clr2.txt", row.names = F, quote = F)

# remove all taxa below 0.2
taxa_over_20pct <- f1y.GBM.clr2 %>% 
  dplyr::filter(prevalence > .2) %>% 
  pull(tax)

# add ABC numbers
f1y_logrelabu_20pct <- f1y_logrelabu %>% dplyr::select(all_of(taxa_over_20pct)) %>% 
  bind_cols(abcno = get_variable(f1y, "ABCNO"), .)

## add phenotypes
abcno_pheno <- read_excel("ABC0232_J45_cox_cross.xlsx", sheet = "Data") %>% 
  dplyr::select(ABCNO, j45_6yr_ever)

## join phenotypes to taxa table & remove ind with NA
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

## get SNP with main effect do GWAS
cmd_1 <- "plink --bfile /home/amelie/nas/PKU_Inter99/merged_kasper/X_merge_1kG_build_QC_clean_postQC" 
cmd_2 <- "--pheno /home/amelie/nas/PKU_Inter99/merged_kasper/phenotypes.txt --ci 0.95"
cmd_3 <- "--logistic --hide-covar --covar /home/amelie/nas/PKU_Inter99/merged_kasper/GWAS/covar_sex.txt --out asthma_main"
system(paste(cmd_1, cmd_2, cmd_3))

assoc <- fread("asthma_main.assoc.logistic") %>% 
  arrange(P)

assoc$P_adjust <- p.adjust(assoc$P, "fdr")
 
snp_main <- assoc %>% 
  dplyr::filter(P_adjust<0.05)

## extract SNPs with p < 0.05 after fdr
write.table(snp_main$SNP, "SNPs_to_extract.txt", quote = F, row.names = F)
system("plink --bfile /home/amelie/nas/PKU_Inter99/merged_kasper/X_merge_1kG_build_QC_clean_postQC --extract SNPs_to_extract.txt --make-bed --out SNPs_main")

setwd("/home/amelie/nas/microbiome/GWAS_main_effects/")

## clump SNPs
system("plink --bfile SNPs_main --clump asthma_main.assoc.logistic --clump-kb 1000 --clump-p1 0.99 --clump-p2 0.99 --clump-r2 0.01 --pheno /home/amelie/nas/PKU_Inter99/merged_kasper/phenotypes.txt --out SNPs_main_clumped") 

## look at error message "missing top variants" in file "SNPs_main_clumped.log"
## is OK, because it has only SNPs with main effects and not all SNPs

## take column from clump file
system("awk '{print $3}' SNPs_main_clumped.clumped > keep_SNPs_clump.txt")

system("plink --bfile /home/amelie/nas/PKU_Inter99/merged_kasper/X_merge_1kG_build_QC_clean_postQC --extract keep_SNPs_clump.txt --make-bed --out SNPs_main_clumped_final")


SNPs_main <- BEDMatrix::BEDMatrix("SNPs_main_clumped_final.bed", simple_names = TRUE)













