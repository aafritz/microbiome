setwd("/home/amelie/microbiome")
rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T, javascript=F)
library(parallel)

data <- fread("f1y_logrelabu_20pct_pheno.txt") 
## 608 individuals: 1 = asthma, 0 = no asthma
## 133 cases already coded in 1
## 475 controls


## if covariates are included they need to be scaled
## phenotype needs to be 0 and 1
## think about scaling input data?

model <- stan_model(file = "/nas/users/amelie/PKU_Inter99/merged_kasper/model/model_anders_nocovar.stan")  
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


