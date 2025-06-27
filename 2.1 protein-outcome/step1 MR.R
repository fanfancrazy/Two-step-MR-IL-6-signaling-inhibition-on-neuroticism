# Load required libraries
library(TwoSampleMR)
library(ggplot2)
library(foreach)
library(data.table)

# Set working directory
setwd("")

# Read list of exposure IDs
iddf = read.table("id.txt", header = TRUE, sep = "\t")
bioid = as.vector(iddf$id)

# Initialize result container
result = data.frame()

# Loop over each exposure ID
foreach(i = bioid, .errorhandling = "pass") %do% {
  
  # Read exposure data
  exposure <- read_exposure_data(
    filename = paste0("cis_decode/", i, "_export.txt"),
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure"
  )
  
  # Read outcome data
  outcome <- read_outcome_data(
    snps = exposure$SNP,
    filename = "ebi-a-GCST006940.txt.gz",
    sep = "\t",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    eaf_col = "eaf.exposure",
    pval_col = "pval.exposure",
    samplesize_col = "samplesize.exposure"
  )
  
  # Harmonize exposure and outcome data
  harmonize <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome, action = 2)
  
  # Filter weak instruments (p-value threshold in outcome)
  harmonize = harmonize[harmonize$pval.outcome > 5e-08,]
  
  # Calculate F-statistic for instrument strength
  harmonize$R2 <- (2 * (harmonize$beta.exposure^2) * harmonize$eaf.exposure * (1 - harmonize$eaf.exposure)) /
    (2 * (harmonize$beta.exposure^2) * harmonize$eaf.exposure * (1 - harmonize$eaf.exposure) +
       2 * harmonize$samplesize.exposure * harmonize$eaf.exposure * (1 - harmonize$eaf.exposure) * harmonize$se.exposure^2)
  
  harmonize$f <- harmonize$R2 * (harmonize$samplesize.exposure - 2) / (1 - harmonize$R2)
  harmonize$meanf <- mean(harmonize$f)
  harmonize <- harmonize[harmonize$f > 10, ]  # Keep only strong instruments
  
  # Run Mendelian Randomization
  mr_result <- mr(harmonize)
  
  # Run Steiger directionality test
  steiger_test <- directionality_test(harmonize)
  
  # Convert effect size to OR (if outcome is binary)
  result_or <- generate_odds_ratios(mr_result)
  
  # Extract p-value
  p_value <- if (length(result_or$pval) >= 3) result_or$pval[3] else result_or$pval[1]
  
  # Save main result
  result = rbind(result, cbind(id = i, pvalue = p_value))
  filename = paste0("result/", i)
  dir.create(filename, showWarnings = FALSE)
  
  # Save harmonized dataset
  write.table(harmonize, file = paste0(filename, "/harmonise.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Save MR result (OR values)
  write.table(result_or[, 5:ncol(result_or)], file = paste0(filename, "/OR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  # Save Steiger test results
  write.table(steiger_test, file = paste0(filename, "/steiger_test.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Pleiotropy test
  pleiotropy = mr_pleiotropy_test(harmonize)
  write.table(pleiotropy, file = paste0(filename, "/pleiotropy.txt"), sep = "\t", quote = FALSE)
  
  # Heterogeneity test
  heterogeneity = mr_heterogeneity(harmonize)
  write.table(heterogeneity, file = paste0(filename, "/heterogeneity.txt"), sep = "\t", quote = FALSE)
  
  # Scatter plot
  p1 <- mr_scatter_plot(mr_result, harmonize)
  ggsave(p1[[1]], file = paste0(filename, "/scatter.pdf"), width = 8, height = 8)
  
  # Single SNP analysis and forest plot
  singlesnp_res <- mr_singlesnp(harmonize)
  singlesnpOR = generate_odds_ratios(singlesnp_res)
  write.table(singlesnpOR, file = paste0(filename, "/singlesnpOR.txt"), row.names = FALSE, sep = "\t", quote = FALSE)
  
  p2 <- mr_forest_plot(singlesnp_res)
  ggsave(p2[[1]], file = paste0(filename, "/forest.pdf"), width = 8, height = 8)
  
  # Leave-one-out analysis
  sen_res <- mr_leaveoneout(harmonize)
  p3 <- mr_leaveoneout_plot(sen_res)
  ggsave(p3[[1]], file = paste0(filename, "/sensitivity-analysis.pdf"), width = 8, height = 8)
  
  # Funnel plot
  p4 <- mr_funnel_plot(singlesnp_res)
  ggsave(p4[[1]], file = paste0(filename, "/funnelplot.pdf"), width = 8, height = 8)
  
  # MR-PRESSO test
  presso = run_mr_presso(harmonize, NbDistribution = 1000)
  capture.output(presso, file = paste0(filename, "/presso.txt"))
}

# Save overall summary result
write.table(result, "result.txt", sep = "\t", quote = FALSE, row.names = FALSE)
