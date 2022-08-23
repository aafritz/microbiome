setwd("/home/amelie/nas/PKU_Inter99/merged_kasper/model/not_stratified/sex_interactions_scaled/")
rm(list=ls())
library(BEDMatrix)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T, javascript=F)
library(parallel)

bed <- BEDMatrix("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_data/clumped_data/PKU_Inter99_clumped.bed", simple_names=TRUE)

fam <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/data/clumped_data/clumped_data/PKU_Inter99_clumped.fam") %>% dplyr::select(FID=V1, IID=V2, PID=V3, MID=V4, SEX=V5, PHENO=V6)

phenotypes <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/phenotypes.txt") %>% dplyr::rename(FID=V1, IID=V2, PHENO=V3)
fam <- left_join(fam, phenotypes, by=c("IID"="IID")) %>% dplyr::select(FID=FID.x, IID, PID, MID, SEX, PHENO=PHENO.y) 

covar <- fread("/home/amelie/nas/PKU_Inter99/merged_kasper/GWAS/covar.txt", header=TRUE)
X <- left_join(fam, covar, by=c("FID"="FID")) %>% dplyr::select(FID, IID.x, PID, MID, SEX, PHENO, C1, C2, C3, C4, C5, C6, C7, C8, C9, C10)
X$SEX <- as.numeric(X$SEX)
covar <- X %>% dplyr::select(C1, C2, C3, C4, C5, C6, C7, C8, C9, C10) %>% as.matrix() %>% scale()

pheno <- X$PHENO 
pheno <- pheno-1 

sex <- X$SEX-1 
snp <- "rs12092254"
model <- stan_model(file = "/nas/users/amelie/PKU_Inter99/merged_kasper/model/model_anders.stan")  
system.time(
  r <- mclapply(colnames(bed), function(snp) {
    
    g_SNP <- bed[,snp] 
    
    X <- cbind(sex, g_SNP) 
    X <- cbind(X, 0)
    
    X <- cbind(X, covar) 
    X <- cbind(X, pheno) %>% na.omit()
    
    X[,"g_SNP"] <- X[,"g_SNP"] - mean(X[,"g_SNP"])
    X[,"sex"] <- X[,"sex"] - mean(X[,"sex"])
    X[,3] <- X[,"g_SNP"] * X[,"sex"]
    
    
    data <- list(N = nrow(X), # number of samples (phenotypes)
                 pheno = X[,"pheno"], # phenotype
                 M = 3, # number of SNPs
                 P = 10, # number of covariates
                 X = X[,1:13]) # predictor matrix
    
    fit <- sampling(
      model,  # Stan program
      data = data,    # named list of data
      chains = 4,             # number of Markov chains
      warmup = 100,          # number of warmup iterations per chain
      iter = 1000,            # total number of iterations per chain
      cores = 4,              # number of cores (could use one per chain)
      refresh = 0             # no progress shown
    )
    
    fit_summary <- summary(fit)
    
    beta_sex <- round(as.numeric(fit_summary$summary[2,"mean"]), 2)
    sd_sex <- round(as.numeric(fit_summary$summary[2,"sd"]), 2)
    n_eff_sex <- round(as.numeric(fit_summary$summary[2,"n_eff"]), 2)
    rhat_sex <- round(as.numeric(fit_summary$summary[2,"Rhat"]), 2)
    
    beta_SNP <- round(as.numeric(fit_summary$summary[3,"mean"]), 2)
    sd_SNP<- round(as.numeric(fit_summary$summary[3,"sd"]), 2)
    n_eff_SNP <- round(as.numeric(fit_summary$summary[3,"n_eff"]), 2)
    rhat_SNP <- round(as.numeric(fit_summary$summary[3,"Rhat"]), 2)
    
    beta_int <- round(as.numeric(fit_summary$summary[4,"mean"]), 2)
    sd_int <- round(as.numeric(fit_summary$summary[4,"sd"]), 2)
    n_eff_int <- round(as.numeric(fit_summary$summary[4,"n_eff"]), 2)
    rhat_int <- round(as.numeric(fit_summary$summary[4,"Rhat"]), 2)
    
    print(c(snp, beta_sex, sd_sex, n_eff_sex, rhat_sex, beta_SNP, sd_SNP, n_eff_SNP, rhat_SNP, 
            beta_int, sd_int, n_eff_int, rhat_int))
    
  }, mc.cores=5)
)

result_df <- as.data.frame(matrix(unlist(r), ncol=13, byrow=TRUE))
colnames(result_df) <- c("snp", "beta_sex", "sd_sex", "n_eff_sex", "rhat_sex", 
                         "beta_SNP", "sd_SNP", "n_eff_SNP", "rhat_SNP", 
                         "beta_int", "sd_int", "n_eff_int", "rhat_int")
write.table(result_df, file = "result_sex_scaled_bayes.txt", quote = FALSE, row.names = FALSE)


