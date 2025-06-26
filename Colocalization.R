# Load required packages
library(dplyr)
library(coloc)
library(data.table)

# Set working directory
setwd("E:/fenland验证/共定位")

# Load pQTL data (exposure dataset)
df1 <- fread("ATF6B_11387_3.txt.gz", header = TRUE)
df1 <- data.frame(
  SNP = df1$rsid,
  chrom = df1$chr,
  pos = df1$pos,
  A1 = toupper(df1$Allele1),
  A2 = toupper(df1$Allele2),
  beta = df1$Effect,
  se = df1$StdErr,
  MAF = df1$Freq1,
  pvalues = df1$Pvalue,
  N = 10708
) %>% filter(!is.na(MAF))

# Load GWAS summary statistics (outcome dataset)
df2 <- fread("ebi-a-GCST006946.txt.gz", header = TRUE)
df2 <- data.frame(
  SNP = df2$SNP,
  chrom = df2$chr.exposure,
  pos = df2$pos.exposure,
  A1 = toupper(df2$effect_allele.exposure),
  A2 = toupper(df2$other_allele.exposure),
  beta = df2$beta.exposure,
  se = df2$se.exposure,
  pvalues = df2$pval.exposure
)

# Define gene region ±250kb
geneChr = 6
geneStart = 32115264
geneEnd = 32128253

df1 <- subset(df1, chrom == geneChr & pos >= geneStart - 250000 & pos <= geneEnd + 250000)
df2 <- subset(df2, chrom == geneChr & pos >= geneStart - 250000 & pos <= geneEnd + 250000)

# Harmonize alleles and merge datasets
dfall <- merge(df1, df2, by = "SNP")
dfall <- dfall[!duplicated(dfall$SNP), ]
dfall <- dfall %>% filter((A1.x == A1.y & A2.x == A2.y) | (A1.x == A2.y & A2.x == A1.y))
dfall <- dfall %>% mutate(beta.y = ifelse(A1.x == A1.y, beta.y, -beta.y))

# Calculate variance
dfall$VAR1 <- dfall$se.x^2
dfall$VAR2 <- dfall$se.y^2
dfall <- dfall[dfall$VAR1 != 0 & dfall$VAR2 != 0, ]

# Format input for coloc
cdf1 <- as.list(data.frame(
  beta = dfall$beta.x,
  varbeta = dfall$VAR1,
  snp = dfall$SNP,
  MAF = dfall$MAF,
  N = dfall$N,
  type = "quant"
))
cdf2 <- as.list(data.frame(
  beta = dfall$beta.y,
  varbeta = dfall$VAR2,
  snp = dfall$SNP,
  type = "cc"
))

# Run coloc analysis
colocresult <- coloc.abf(cdf1, cdf2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
coloc_summary <- colocresult$summary
coloc_results <- colocresult$results

# Save results
write.csv(coloc_results, "coloc_results.csv", row.names = FALSE)
write.csv(as.data.frame(t(coloc_summary)), "coloc_summary.csv", row.names = FALSE)

# Add -log10(P) for plotting
dfall$P1 <- 2 * pnorm(-abs(dfall$beta.x / sqrt(dfall$VAR1)))
dfall$P2 <- 2 * pnorm(-abs(dfall$beta.y / sqrt(dfall$VAR2)))
dfall$logP1 <- -log10(dfall$P1)
dfall$logP2 <- -log10(dfall$P2)

# Merge PP4 results for plotting
dfplot <- merge(dfall, coloc_results[, c("snp", "SNP.PP.H4")], by.x = "SNP", by.y = "snp", all.x = TRUE)

# -----------------------------
# Optional: Locuscompare plot
# -----------------------------
library(gwasglue)
library(locuscomparer)
library(openxlsx)
library(tidyverse)

pdf("locuscompare_plot3.pdf")

gwas_fn <- df2[, c("SNP", "pvalues")] %>%
  rename(rsid = SNP, pval = pvalues)
pqtl_fn <- df1[, c("SNP", "pvalues")] %>%
  rename(rsid = SNP, pval = pvalues)

print(locuscompare(
  in_fn1 = gwas_fn,
  in_fn2 = pqtl_fn,
  title1 = "GWAS",
  title2 = "pQTL"
))
dev.off()
