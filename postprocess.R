setwd("/home/amelie/nas/microbiome")
rm(list=ls())

res <- fread("res_bayes.txt") %>% 
  dplyr::select(taxa, beta_int, sd_int) %>% 
  mutate(prob=pnorm(0.1, abs(beta_int), sd_int, lower.tail=FALSE), OR=exp(beta_int)) %>% 
  arrange(desc(prob))
