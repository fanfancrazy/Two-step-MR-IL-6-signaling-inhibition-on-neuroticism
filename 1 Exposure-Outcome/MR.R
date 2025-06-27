# === Load necessary libraries ===
library(TwoSampleMR)

# === Set parameters ===
exposureName <- "CRP Level || IL6R"   # 用于注释
diseaseName <- "CRP"                  # 用于注释
outcomeID <- "ebi-a-GCST006940"       # 在线结局ID
setwd("\\IL6-Neurociticism\\IL6R")      # 修改为你的实际路径

# === 从 Supplementary 文件读取暴露数据 ===
exposure_dat <- read_exposure_data(
  filename = "Supplementary Table2",
  #filename = "Supplementary Table4",
  sep = ",",
  snp_col = "SNP",
  beta_col = "beta.exposure",
  se_col = "se.exposure",
  effect_allele_col = "effect_allele.exposure",
  other_allele_col = "other_allele.exposure",
  eaf_col = "eaf.exposure",  # 如无可删去此行
  pval_col = "pval.exposure",
  samplesize_col = "samplesize.exposure"  # 如无可删去此行
)

# === 从在线数据库提取结局数据 ===
outcome_dat <- extract_outcome_data(
  snps = exposure_dat$SNP,
  outcomes = outcomeID
)

# === Harmonize ===
exposure_dat$exposure <- exposureName
outcome_dat$outcome <- diseaseName
dat <- harmonise_data(exposure_dat, outcome_dat)
dat <- dat[dat$pval.outcome > 5e-08,]

# === Keep SNPs suitable for MR analysis and export ===
outTab <- subset(dat, mr_keep == TRUE)
write.csv(outTab, file = "table.SNP.csv", row.names = FALSE)

# === Perform MR-PRESSO to detect potential outliers ===
presso <- run_mr_presso(dat)
write.csv(presso[[1]]$`MR-PRESSO results`$`Global Test`, file = "table.MR-PRESSO_Global.csv")
write.csv(presso[[1]]$`MR-PRESSO results`$`Outlier Test`, file = "table.MR-PRESSO_Outlier.csv")

# === Perform Mendelian Randomization analysis ===
mrResult <- mr(dat)
mrTab <- generate_odds_ratios(mrResult)
write.csv(mrTab, file = "table.MRresult.csv", row.names = FALSE)

# === Heterogeneity test ===
heterTab <- mr_heterogeneity(dat)
write.csv(heterTab, file = "table.heterogeneity.csv", row.names = FALSE)

# === Pleiotropy test ===
pleioTab <- mr_pleiotropy_test(dat)
write.csv(pleioTab, file = "table.pleiotropy.csv", row.names = FALSE)

# === Plot: Scatter plot of MR estimates ===
pdf(file = "pic.scatter_plot.pdf", width = 7.5, height = 7)
mr_scatter_plot(mrResult, dat)
dev.off()

# === Plot: Forest plot of individual SNP effects ===
res_single <- mr_singlesnp(dat)
pdf(file = "pic.forest.pdf", width = 7, height = 5.5)
mr_forest_plot(res_single)
dev.off()

# === Plot: Funnel plot for detecting asymmetry ===
pdf(file = "pic.funnel_plot.pdf", width = 7, height = 6.5)
mr_funnel_plot(singlesnp_results = res_single)
dev.off()

# === Plot: Leave-one-out sensitivity analysis ===
pdf(file = "pic.leaveoneout.pdf", width = 7, height = 5.5)
mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
dev.off()


