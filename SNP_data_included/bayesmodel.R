setwd("/home/amelie/nas/microbiome/SNP_data_included/")
rm(list=ls())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = T, javascript=F)
library(parallel)

data <- fread("f1y_g.GBM.clr2_20pct_t_abcno_pheno_sex_final_gt.txt") 
# 1 = asthma, 0 = no asthma


## if covariates are included they need to be scaled
## phenotype needs to be 0 and 1
## think about scaling input data?
sex <- data$SEX
pheno <- data$pheno
i <- 4
model <- stan_model(file = "model_nocovar_taxa_sex_snp.stan")  
system.time(
  r <- mclapply(4:ncol(data), function(i) {
#  r <- mclapply(4:10, function(i) {    
    taxa <- data[,..i] 
    taxa_name <- colnames(taxa)
    
    X <- cbind(sex, taxa) 
    X <- cbind(X, 0)
    colnames(X)[2:3] <- c("taxa", "int")
    X <- cbind(X, data[,62:ncol(data)])
    
    X <- cbind(X, pheno) %>% 
      na.omit() 
    
    X[,"sex"] <- X[,"sex"] - mean(X$sex)
    X[,"int"] <- X[,"taxa"] * X[,"sex"]
    
    data <- list(N = nrow(X), # number of samples (phenotypes)
                 pheno = X$pheno, # phenotype
                 M = 55, # sex, taxa, interaction, 52 SNPs
                 #P = 10, # number of covariates
                 X = X[,1:55]) # predictor matrix
    
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
    
    beta_taxa <- round(as.numeric(fit_summary$summary[3,"mean"]), 2)
    sd_taxa<- round(as.numeric(fit_summary$summary[3,"sd"]), 2)
    n_eff_taxa <- round(as.numeric(fit_summary$summary[3,"n_eff"]), 2)
    rhat_taxa <- round(as.numeric(fit_summary$summary[3,"Rhat"]), 2)
    
    beta_int <- round(as.numeric(fit_summary$summary[4,"mean"]), 2)
    sd_int <- round(as.numeric(fit_summary$summary[4,"sd"]), 2)
    n_eff_int <- round(as.numeric(fit_summary$summary[4,"n_eff"]), 2)
    rhat_int <- round(as.numeric(fit_summary$summary[4,"Rhat"]), 2)
    
    print(c(taxa_name, beta_sex, sd_sex, n_eff_sex, rhat_sex, beta_taxa, sd_taxa, n_eff_taxa, rhat_taxa, 
            beta_int, sd_int, n_eff_int, rhat_int))
    
  }, mc.cores=20)
)

result_df <- as.data.frame(matrix(unlist(r), ncol=13, byrow=TRUE))
colnames(result_df) <- c("taxa", "beta_sex", "sd_sex", "n_eff_sex", "rhat_sex", 
                         "beta_taxa", "sd_taxa", "n_eff_taxa", "rhat_taxa", 
                         "beta_int", "sd_int", "n_eff_int", "rhat_int")
write.table(result_df, file = "res_bayes.txt", quote = FALSE, row.names = FALSE)


