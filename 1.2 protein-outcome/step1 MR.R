# Load necessary libraries
library(TwoSampleMR)
library(ggplot2)
library(foreach)

# Set working directory (change this to your local path if needed)
setwd(":\\Colocalization\\IL6_protein_project")

# Load outcome trait IDs to be looped over
iddf <- read.table("id.txt", header = TRUE, sep = "\t")
bioid <- as.vector(iddf$id)

# Initialize result container
result <- data.frame()

# Loop over each outcome trait ID
foreach(i = bioid, .errorhandling = "pass") %do% {
  
  # Load exposure data (e.g., IL6R instruments)
  expo_rt <- read_exposure_data(
    filename = "IL6R_instruments_negated.csv",  # e.g., beta flipped
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
  
  # Load corresponding outcome data
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
  
  # Harmonize exposure and outcome datasets
  harm_rt <- harmonise_data(exposure_dat = expo_rt, outcome_dat = outc_rt, action = 2)
  
  # Optional: filter out strongly associated outcomes (P > 5e-8)
  harm_rt <- harm_rt[harm_rt$pval.outcome > 5e-8, ]
  
  # Calculate F-statistics
  harm_rt$R2 <- (2 * harm_rt$beta.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure)) /
    (2 * harm_rt$beta.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure) +
       2 * harm_rt$samplesize.exposure * harm_rt$se.exposure^2 * harm_rt$eaf.exposure * (1 - harm_rt$eaf.exposure))
  
  harm_rt$f <- harm_rt$R2 * (harm_rt$samplesize.exposure - 2) / (1 - harm_rt$R2)
  harm_rt$meanf <- mean(harm_rt$f, na.rm = TRUE)
  
  # Filter weak instruments (F > 15)
  harm_rt <- harm_rt[harm_rt$f > 15, ]
  
  # Perform MR analysis
  mr_result <- mr(harm_rt)
  result_or <- generate_odds_ratios(mr_result)
  
  # If results are valid (p-value not NA), proceed to save
  if (!is.na(mr_result$pval[3]) && mr_result$pval[3] < 1) {
    
    result <- rbind(result, cbind(id = i, pvalue = result_or$pval[3]))
    
    # Create a result folder for this outcome
    filename <- paste0("result/", i)
    dir.create(filename, recursive = TRUE)
    
    # Save harmonized data and MR results
    write.table(harm_rt, file = paste0(filename, "/harmonise.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Steiger directionality test
    C <- directionality_test(harm_rt)
    write.table(C, file = paste0(filename, "/steiger_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    # Pleiotropy and heterogeneity tests
    pleiotropy <- mr_pleiotropy_test(harm_rt)
    heterogeneity <- mr_heterogeneity(harm_rt)
    write.table(pleiotropy, file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)
    write.table(heterogeneity, file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)
    
    # Plot MR scatter plot
    p1 <- mr_scatter_plot(mr_result, harm_rt)
    ggsave(p1[[1]], file = paste0(filename, "/scatter.pdf"), width = 8, height = 8)
    
    # Forest plot and single SNP analysis
    singlesnp_res <- mr_singlesnp(harm_rt)
    singlesnpOR <- generate_odds_ratios(singlesnp_res)
    write.table(singlesnpOR, file = paste0(filename, "/singlesnpOR.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
    
    p2 <- mr_forest_plot(singlesnp_res)
    ggsave(p2[[1]], file = paste0(filename, "/forest.pdf"), width = 8, height = 8)
    
    # Leave-one-out sensitivity analysis
    sen_res <- mr_leaveoneout(harm_rt)
    p3 <- mr_leaveoneout_plot(sen_res)
    ggsave(p3[[1]], file = paste0(filename, "/sensitivity-analysis.pdf"), width = 8, height = 8)
    
    # Funnel plot
    p4 <- mr_funnel_plot(singlesnp_res)
    ggsave(p4[[1]], file = paste0(filename, "/funnelplot.pdf"), width = 8, height = 8)
    
    # MR-PRESSO analysis
    presso <- run_mr_presso(harm_rt, NbDistribution = 1000)
    capture.output(presso, file = paste0(filename, "/presso.txt"))
  }
}

# Save summary p-values across all outcomes
write.table(result, "r_result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
