# Load necessary library
library(TwoSampleMR)

# === Set basic parameters ===
exposureID <- "ieu-b-35"              # ID for exposure dataset (CRP level)
outcomeID <- "ebi-a-GCST006940"       # ID for outcome dataset
geneChr <- 1                          # Chromosome number where IL6R is located
geneStart <- 154377819                # Start position of IL6R gene
geneEnd <- 154441926                  # End position of IL6R gene

exposureName <- "CRP Level || IL6R"   # Label for exposure in plots/tables
diseaseName <- "CRP"                  # Label for outcome
setwd("/IL6-Neurociticism/IL6R")    # Set working directory (modify as needed)

# === Extract and clump exposure instruments ===
exposure_dat <- extract_instruments(outcome = exposureID, clump = FALSE)
exposure_dat <- clump_data(exposure_dat, clump_kb = 100, clump_r2 = 0.3)

# === Filter SNPs within ±100kb of the IL6R gene region and EAF > 0.01 ===
geneData <- subset(exposure_dat, chr.exposure == geneChr)
geneData <- subset(geneData, pos.exposure >= (geneStart - 1e5) & pos.exposure <= (geneEnd + 1e5))
geneData <- subset(geneData, eaf.exposure > 0.01)

# === Calculate F-statistics to assess instrument strength ===
geneData$R2 <- (2 * geneData$beta.exposure^2 * geneData$eaf.exposure * (1 - geneData$eaf.exposure)) /
  (2 * geneData$beta.exposure^2 * geneData$eaf.exposure * (1 - geneData$eaf.exposure) +
     2 * geneData$se.exposure^2 * geneData$samplesize.exposure * geneData$eaf.exposure * (1 - geneData$eaf.exposure))

geneData$F <- geneData$R2 * (geneData$samplesize.exposure - 2) / (1 - geneData$R2)

# === Keep only strong instruments (F > 10) and save to file ===
filtered_exposure <- subset(geneData, F > 10)
write.csv(filtered_exposure, file = "IL6R_instruments.csv", row.names = FALSE)

# === Extract outcome data for selected SNPs ===
outcome_dat <- extract_outcome_data(snps = geneData$SNP, outcomes = outcomeID)

# === Harmonize exposure and outcome data ===
geneData$exposure <- exposureName
outcome_dat$outcome <- diseaseName
dat <- harmonise_data(geneData, outcome_dat)

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
