# Load required libraries
library(TwoSampleMR)
library(ggplot2)
library(foreach)

# Set working directory
setwd("")

# Load list of outcome IDs
iddf <- read.table("id.txt", header = TRUE, sep = "\t")
bioid <- as.vector(iddf$id)
result <- data.frame()

# Loop over each outcome ID
foreach(i = bioid, .errorhandling = "pass") %do% {
  
  # Load exposure data (IL6R instruments)
  exposure <- read_exposure_data(
    filename = "IL6R_instruments.csv",
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
  outcome <- read_outcome_data(
    snps = exposure$SNP,
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
  harmonise <- harmonise_data(exposure_dat = exposure, outcome = outcome, action = 2)
  
  # Filter SNPs with outcome p > 5e-8
  harmonise <- harmonise[harmonise$pval.outcome > 5e-8, ]
  
  # Calculate F-statistics for instrument strength
  harmonise$R2 <- (2 * harmonise$beta.exposure^2 * harmonise$eaf.exposure * (1 - harmonise$eaf.exposure)) /
    (2 * harmonise$beta.exposure^2 * harmonise$eaf.exposure * (1 - harmonise$eaf.exposure) +
       2 * harmonise$samplesize.exposure * harmonise$se.exposure^2 * harmonise$eaf.exposure * (1 - harmonise$eaf.exposure))
  harmonise$f <- harmonise$R2 * (harmonise$samplesize.exposure - 2) / (1 - harmonise$R2)
  harmonise$meanf <- mean(harmonise$f, na.rm = TRUE)
  
  # Keep SNPs with F > 15
  harmonise <- harmonise[harmonise$f > 15, ]
  
  # Run MR analysis
  mr_result <- mr(harmonise)
  result_or <- generate_odds_ratios(mr_result)
  
  if (!is.na(mr_result$pval[3]) && mr_result$pval[3] < 1) {
    result <- rbind(result, cbind(id = i, pvalue = result_or$pval[3]))
    filename <- paste0("result/", i)
    dir.create(filename, showWarnings = FALSE, recursive = TRUE)
    
    # Save harmonised data and OR result
    write.table(harmonise, file = paste0(filename, "/harmonise.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Sensitivity analyses
    write.table(directionality_test(harmonise), file = paste0(filename, "/steiger_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(mr_pleiotropy_test(harmonise), file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)
    write.table(mr_heterogeneity(harmonise), file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)
    
    # Plots
    ggsave(mr_scatter_plot(mr_result, harmonise)[[1]], file = paste0(filename, "/scatter.pdf"), width = 8, height = 8)
    singlesnp_res <- mr_singlesnp(harmonise)
    write.table(generate_odds_ratios(singlesnp_res), file = paste0(filename, "/singlesnpOR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    ggsave(mr_forest_plot(singlesnp_res)[[1]], file = paste0(filename, "/forest.pdf"), width = 8, height = 8)
    sen_res <- mr_leaveoneout(harmonise)
    ggsave(mr_leaveoneout_plot(sen_res)[[1]], file = paste0(filename, "/sensitivity-analysis.pdf"), width = 8, height = 8)
    ggsave(mr_funnel_plot(singlesnp_res)[[1]], file = paste0(filename, "/funnelplot.pdf"), width = 8, height = 8)
    
    # MR-PRESSO analysis
    presso <- run_mr_presso(harmonise, NbDistribution = 1000)
    capture.output(presso, file = paste0(filename, "/presso.txt"))
  }
}

# Save summary p-value table
write.table(result, "r_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
