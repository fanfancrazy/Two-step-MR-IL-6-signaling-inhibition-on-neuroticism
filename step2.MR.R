# ================================================
# Batch Mendelian Randomization for IL6R Instruments
# ================================================
# This script performs two-sample MR across multiple outcomes
# using IL6R instruments and saves harmonised data, MR results,
# and various sensitivity test outputs for each outcome.
# ================================================

library(TwoSampleMR)
library(ggplot2)
library(foreach)

# Set working directory
setwd("E:/Colocalization/IL6_protein_project")

# Load list of outcome IDs
iddf <- read.table("id.txt", header = TRUE, sep = "\t")
bioid <- as.vector(iddf$id)
result <- data.frame()

# Loop over each outcome ID
foreach(i = bioid, .errorhandling = "pass") %do% {
  
  # Load exposure data (IL6R instruments)
  expo_rt <- read_exposure_data(
    filename = "IL6R_instruments_取负号.csv",
    sep = ",",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure"
  )
  
  # Load outcome data
  outc_rt <- read_outcome_data(
    snps = expo_rt$SNP,
    filename = paste0("multi_outcome/", i, ".txt.gz"),
    sep = "\t",
    snp_col = "rsids",
    beta_col = "Beta",
    se_col = "SE",
    effect_allele_col = "effectAllele",
    other_allele_col = "otherAllele",
    eaf_col = "ImpMAF",
    pval_col = "Pval",
    samplesize_col = "N"
  )
  
  # Harmonise datasets
  harm_rt <- harmonise_data(exposure_dat = expo_rt, outcome_dat = outc_rt, action = 2)
  
  # Filter SNPs with outcome p > 5e-8
  harm_rt <- harm_rt[harm_rt$pval.outcome > 5e-8, ]
  
  # Calculate F-statistics for instrument strength
  harm_rt$R2 <- (2 * harm_rt$beta.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * harm_rt$beta.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure))
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf <- mean(harm_rt$f, na.rm = TRUE)
  
  # Keep SNPs with F > 15
  harm_rt <- harm_rt[harm_rt$f > 15, ]
  
  # Run MR analysis
  mr_result <- mr(harm_rt)
  result_or <- generate_odds_ratios(mr_result)
  
  if (!is.na(mr_result$pval[3]) && mr_result$pval[3] < 1) {
    result <- rbind(result, cbind(id = i, pvalue = result_or$pval[3]))
    filename <- paste0("result/", i)
    dir.create(filename, showWarnings = FALSE, recursive = TRUE)
    
    # Save harmonised data and OR result
    write.table(harm_rt, file = paste0(filename, "/harmonise.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Sensitivity analyses
    write.table(directionality_test(harm_rt), file = paste0(filename, "/steiger_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(mr_pleiotropy_test(harm_rt), file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)
    write.table(mr_heterogeneity(harm_rt), file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)
    
    # Plots
    ggsave(mr_scatter_plot(mr_result, harm_rt)[[1]], file = paste0(filename, "/scatter.pdf"), width = 8, height = 8)
    singlesnp_res <- mr_singlesnp(harm_rt)
    write.table(generate_odds_ratios(singlesnp_res), file = paste0(filename, "/singlesnpOR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(mr_forest_plot(singlesnp_res)[[1]], file = paste0(filename, "/forest.pdf"), width = 8, height = 8)
    sen_res <- mr_leaveoneout(harm_rt)
    ggsave(mr_leaveoneout_plot(sen_res)[[1]], file = paste0(filename, "/sensitivity-analysis.pdf"), width = 8, height = 8)
    ggsave(mr_funnel_plot(singlesnp_res)[[1]], file = paste0(filename, "/funnelplot.pdf"), width = 8, height = 8)
    
    # MR-PRESSO analysis
    presso <- run_mr_presso(harm_rt, NbDistribution = 1000)
    capture.output(presso, file = paste0(filename, "/presso.txt"))
  }
}

# Save summary p-value table
write.table(result, "r_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
