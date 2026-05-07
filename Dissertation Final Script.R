# 1. Load Packages 

library(SummarizedExperiment)
library(HDF5Array)
library(DelayedArray)
library(pheatmap)
library(ggplot2)
library(survival)
library(survminer)
library(DelayedMatrixStats)

# 2. Load Saved HDF5 Objects 

tcga_common_h5 <- loadHDF5SummarizedExperiment("tcga_brca_common_hdf5")
metabric_common_h5 <- loadHDF5SummarizedExperiment("metabric_common_hdf5")

# METABRIC assay may be unnamed
if (is.null(assayNames(metabric_common_h5))) {
  assayNames(metabric_common_h5) <- "data"
}

assay_tcga <- "unstranded"
assay_meta <- assayNames(metabric_common_h5)[1]


# 3. PAM50 Genes 

pam50 <- c(
  "BAG1","BCL2","BIRC5","BLGN","CCNB1","CDC20","CDC6","CDH3","CENPF","CEP55",
  "CXXC5","EGFR","ERBB2","ESR1","EXO1","FGFR4","FOXA1","GRP160","GRB7","KIF2C",
  "KRT14","KRT17","KRT5","KRT8","KRT18","KRT19","MELK","MIA","MK167","MLPH",
  "MYBL2","MYC","NAT1","ORC6L","PGR","PHGDH","PTTG1","RRM2","SLC39A6","TMEM45B"
  ,"TYMS","UBE2C","BCL2L2","KRT7","KIF20A","MMP11","CSTA","GSTM1","PAMR1","KRT6B"
)

pam50_overlap <- intersect(toupper(pam50), rownames(tcga_common_h5))
pam50_overlap <- intersect(pam50_overlap, rownames(metabric_common_h5))

# 4. PAM50 Signature Scoring

score_signature <- function(se, assay_name, genes) {
  X <- as.matrix(assay(se[genes, ], i= assay_name, withDimnames = TRUE))
  X <- log2(X + 1)
  Xz <- t(scale(t(X)))     # z-score per gene
  colMeans(Xz, na.rm=TRUE) # mean z across genes = signature score
}

pam50_score_tcga  <- score_signature(tcga_common_h5, assay_tcga, pam50_overlap)
pam50_score_meta  <- score_signature(metabric_common_h5, assay_meta, pam50_overlap)

saveRDS(pam50_score_tcga, "pam50_score_tcga.rds")
saveRDS(pam50_score_meta, "pam50_score_metabric.rds")

# 5. PAM50 Distribution Plots 

boxplot(pam50_score_tcga, pam50_score_meta,
        names = c("TCGA-BRCA","METABRIC"),
        ylab = "PAM50 (mean z-score)",
        main = "PAM50 signature score distribution")

# Density Plot 
d1 <- density(pam50_score_tcga)
d2 <- density(pam50_score_meta)

xlim <- range(c(d1$x, d2$x))
ylim <- range(c(d1$y, d2$y))

plot(d1, xlim = xlim, ylim = ylim,
     main = "PAM50 score distribution",
     xlab = "Score",
     col = "steelblue", lwd = 2)
lines(d2, col = "firebrick", lwd = 2)
legend("topleft", legend = c("TCGA-BRCA", "METABRIC"),
       col = c("steelblue","firebrick"), lwd = 2, bty = "n")


# 6. PAM50 Heatmaps 

mat_tcga_z <- log2(as.matrix(assay(tcga_common_h5, assay_tcga)[pam50_overlap, ]) + 1)
mat_meta_z <- log2(as.matrix(assay(metabric_common_h5, assay_meta)[pam50_overlap, ]) + 1)

par(mar = c(8,10, 8, 8))
pheatmap(mat_tcga_z,
         main = "TCGA PAM50 Gene Expression Heatmap (z-scored)",
         show_colnames = FALSE,
         cellheight = 6,
         fontsize_row = 6)

par(mar = c(8, 10, 8 , 8))
pheatmap(mat_meta_z,
         main = "METABRIC PAM50 Gene Expression Heatmap (z-scored)",
         show_colnames = FALSE,
         cellheight = 6,
         fontsize_row = 6)


# 7. PAM50 Concordance Analysis 

pam50_overlap <- intersect(toupper(pam50), rownames(tcga_common_h5))
pam50_overlap <- intersect(pam50_overlap, rownames(metabric_common_h5))

tcga_mat <- log2(as.matrix(assay(tcga_common_h5)[pam50_overlap, ]) + 1)
meta_mat <- log2(as.matrix(assay(metabric_common_h5)[pam50_overlap, ]) + 1)

tcga_gene_means <- rowMeans(tcga_mat)
meta_gene_means <- rowMeans(meta_mat)

pam50_corr <- cor(tcga_gene_means, meta_gene_means, method="spearman")
pam50_corr

op <- par(mar = c(5, 5, 4, 12), xpd = FALSE)  # xpd FALSE keeps lines clipped

plot(tcga_gene_means, meta_gene_means,
     xlab = "Mean PAM50 expression (TCGA)",
     ylab = "Mean PAM50 expression (METABRIC)",
     main = paste("PAM50 concordance: ρ =", round(pam50_corr, 2)),
     pch = 16)

# identity and regression line
abline(0, 1, col = "red", lty = 2, lwd = 2)
fit <- lm(meta_gene_means ~ tcga_gene_means)
abline(fit, col = "blue", lwd = 2)

par(xpd = NA)

usr <- par("usr")  # plot region limits: x1, x2, y1, y2
x_out <- usr[2] + 0.03 * (usr[2] - usr[1])  # just outside the plot
y_top <- usr[4]

legend(x = x_out, y = y_top,
       legend = c("Identity line", "Regression line"),
       col = c("red", "blue"),
       lty = c(2, 1),
       lwd = 2,
       bty = "n",
       xjust = 0, yjust = 1)

par(op)


# 8. Basal/ Luminal Module Scores 

luminal_genes <- c("ESR1","PGR","FOXA1","BCL2")
basal_genes   <- c("KRT5","KRT14","KRT17","EGFR")

luminal_genes <- intersect(luminal_genes, pam50_overlap)
basal_genes   <- intersect(basal_genes, pam50_overlap)

lum_tcga <- colMeans(t(scale(t(tcga_mat[luminal_genes, ]))))
bas_tcga <- colMeans(t(scale(t(tcga_mat[basal_genes, ]))))

lum_meta <- colMeans(t(scale(t(meta_mat[luminal_genes, ]))))
bas_meta <- colMeans(t(scale(t(meta_mat[basal_genes, ]))))

boxplot(lum_tcga, lum_meta,
        names=c("TCGA","METABRIC"),
        main="Luminal module score comparison",
        xlab = "Cohort",
        ylab = "Module Score")

text(
  x = 1:2,
  y = boxplot(lum_tcga, lum_meta, plot = FALSE)$stats[3, ],
  labels = round(boxplot(lum_tcga, lum_meta, plot = FALSE)$stats[3, ], 2),
  pos = 3,
  cex = 0.8
)

boxplot(bas_tcga, bas_meta,
        names=c("TCGA","METABRIC"),
        main="Basal module score comparison",
        xlab = "Cohort",
        ylab = "Module Score")

text(
  x = 1:2,
  y = boxplot(bas_tcga, bas_meta, plot = FALSE)$stats[3, ],
  labels = round(boxplot(bas_tcga, bas_meta, plot = FALSE)$stats[3, ], 2),
  pos = 3,
  cex = 0.8 
)


# 9. PCA

pca_tcga <- prcomp(t(tcga_mat), scale.=TRUE)
pca_meta <- prcomp(t(meta_mat), scale.=TRUE)

plot(pca_tcga$x[,1], pca_tcga$x[,2],
     main="TCGA PCA (PAM50 genes)",
     xlab="PC1", ylab="PC2")

plot(
  pca_meta$x[,1], pca_meta$x[,2],
  main="METABRIC PCA (PAM50 genes)",
  xlab="PC1", ylab="PC2"
)


# 10. PAM50 Outlier Table 

delta <- meta_gene_means - tcga_gene_means
abs_delta <- abs(delta)

top_n <- 10   # number of outliers 
top_genes <- names(sort(abs_delta, decreasing = TRUE))[1:top_n]

outlier_table <- data.frame(
  Gene = top_genes,
  TCGA_mean = tcga_gene_means[top_genes],
  METABRIC_mean = meta_gene_means[top_genes],
  Difference_META_minus_TCGA = delta[top_genes],
  Abs_difference = abs_delta[top_genes],
  Direction = ifelse(delta[top_genes] > 0, "Higher in METABRIC", "Higher in TCGA"),
  row.names = NULL
)

outlier_table <- outlier_table[order(outlier_table$Abs_difference, decreasing = TRUE), ]
outlier_table[, c("TCGA_mean","METABRIC_mean","Difference_META_minus_TCGA","Abs_difference")] <-
  round(outlier_table[, c("TCGA_mean","METABRIC_mean","Difference_META_minus_TCGA","Abs_difference")], 3)

outlier_table

write.csv(outlier_table, "PAM50_outlier_genes.csv", row.names = FALSE)


# 11. PAM50 Survival Analysis: TCGA 

cd <- as.data.frame(colData(tcga_common_h5))

surv_tcga_pam50 <- cd
surv_tcga_pam50$time <- ifelse(
  !is.na(surv_tcga_pam50$days_to_death),
  surv_tcga_pam50$days_to_death,
  surv_tcga_pam50$days_to_last_follow_up 
)

surv_tcga_pam50$event <- ifelse(surv_tcga_pam50$vital_status == "Dead", 1, 0)
surv_tcga_pam50$time_years <- surv_tcga_pam50$time / 365.25
surv_tcga_pam50$pam50_score <- pam50_score_tcga[rownames(surv_tcga_pam50)]

surv_tcga_pam50 <- surv_tcga_pam50[
  !is.na(surv_tcga_pam50$time_years) &
  !is.na(surv_tcga_pam50$event) &
  !is.na(surv_tcga_pam50$pam50_score),
]

median_tcga <- median(surv_tcga_pam50$pam50_score, na.rm = TRUE)
surv_tcga_pam50$group <- ifelse(
  surv_tcga_pam50$pam50_score >= median_tcga,
  "High PAM50" ,
  "Low PAM50"
)

fit_tcga_pam50 <- survfit(Surv(time_years, event) ~ group, data = surv_tcga_pam50)
ggsurvplot(
  fit_tcga_pam50,
  data = surv_tcga_pam50,
  palette = c("#0072B2", "#CC79A7"),
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  title = "Kaplan-Meier Overall Survival Stratified by PAM50 Score (TCGA-BRCA)",
  legend.title = "PAM50 Score",
  legend.labs = c("High", "Low"),
  xlab = "Time (years)",
  ylab = "Overall survival probability"
)

cox_tcga_pam50 <- coxph(Surv(time_years, event) ~ pam50_score, data = surv_tcga_pam50)
summary(cox_tcga_pam50)


# 12. PAM50 Survival Analysis: METABRIC

meta_clin <- read.delim("brca_metabric_clinical_data.tsv", comment.char = "#")
status_clean <- trimws(meta_clin$Overall.Survival.Status)

surv_meta_pam50 <- meta_clin 
surv_meta_pam50$time_years <- as.numeric(surv_meta_pam50$Overall.Survival..Months.) / 12 

surv_meta_pam50$event <- ifelse(
  grepl("DECEASED|DEAD", status_clean, ignore.case = TRUE),
  1,
  ifelse(
    grepl("LIVING|ALIVE", status_clean, ignore.case = TRUE),
    0,
    NA
  )
)

surv_meta_pam50$pam50_score <- pam50_score_meta[surv_meta_pam50$Patient.ID]

surv_meta_pam50 <- surv_meta_pam50[
  !is.na(surv_meta_pam50$time_years) &
  !is.na(surv_meta_pam50$event) &
  !is.na(surv_meta_pam50$pam50_score),
]

median_meta <- median(surv_meta_pam50$pam50_score, na.rm = TRUE)
surv_meta_pam50$group <- ifelse(
  surv_meta_pam50$pam50_score >= median_meta,
  "High PAM50",
  "Low PAM50"
)

table(surv_meta_pam50$group)

fit_meta_pam50 <- survfit(Surv(time_years, event) ~ group, data = surv_meta_pam50)

ggsurvplot(
  fit_meta_pam50,
  data = surv_meta_pam50,
  pval = TRUE,
  palette = c("#0072B2", "#CC79A7"),
  risk.table = TRUE,
  conf.int = FALSE,
  title = "Kaplan-Meier Overall Survival Stratified by PAM50 Score (METABRIC)",
  xlab = "Time (years)",
  ylab = "Overall survival probability",
  legend.title = "PAM50 score",
  legend.labs = c("High", "Low")
)

cox_meta_pam50 <- coxph(Surv(time_years, event) ~ pam50_score, data = surv_meta_pam50)


# 13. PAM50 Cox Results Table 

extract_cox_summary <- function(fit, cohort) {
  s <- summary(fit)
  
  data.frame(
    Cohort = cohort,
    HR = s$conf.int[1, "exp(coef)"],
    CI_lower = s$conf.int[1, "lower .95"],
    CI_upper = s$conf.int[1, "upper .95"],
    p_value = s$coefficients[1, "Pr(>|z|)"],
    C_index = s$concordance[1]
  )
}

table_pam50 <- rbind(
  extract_cox_summary(cox_tcga_pam50, "TCGA-BRCA"),
  extract_cox_summary(cox_meta_pam50, "METABRIC")
)

table_pam50[, -1] <- round(table_pam50[, -1], 3)

write.csv(table_pam50, "pam50_cox_results.csv", row.names = FALSE)


# Objective 3: Build survival tables

# TCGA
cd_tcga <- as.data.frame(colData(tcga_common_h5))

surv_tcga <- data.frame(
  sample_id = rownames(cd_tcga),
  time_years = ifelse(
    !is.na(cd_tcga$days_to_death),
    cd_tcga$days_to_death,
    cd_tcga$days_to_last_follow_up
  ) / 365.25,
  event = ifelse(cd_tcga$vital_status == "Dead", 1, 0)
)

surv_tcga <- surv_tcga[
  !is.na(surv_tcga$time_years) &
    !is.na(surv_tcga$event),
]

# METABRIC
tmp_meta <- meta_clin
status_clean <- trimws(tmp_meta$Overall.Survival.Status)

tmp_meta$time_years <- as.numeric(tmp_meta$Overall.Survival..Months.) / 12
tmp_meta$event <- ifelse(
  grepl("DECEASED|DEAD", status_clean, ignore.case = TRUE),
  1,
  ifelse(
    grepl("LIVING|ALIVE", status_clean, ignore.case = TRUE),
    0,
    NA
  )
)

id_col <- if (
  sum(tmp_meta$Sample.ID %in% colnames(metabric_common_h5)) >=
  sum(tmp_meta$Patient.ID %in% colnames(metabric_common_h5))
) "Sample.ID" else "Patient.ID"

surv_meta <- data.frame(
  sample_id = tmp_meta[[id_col]],
  time_years = tmp_meta$time_years,
  event = tmp_meta$event
)

surv_meta <- surv_meta[
  !is.na(surv_meta$sample_id) &
    !is.na(surv_meta$time_years) &
    !is.na(surv_meta$event),
]


# 15. Objective 3: Top Variable Genes

common_genes <- intersect(rownames(tcga_common_h5), rownames(metabric_common_h5))

tcga_log <- log2(assay(tcga_common_h5, i = assay_tcga)[common_genes, ] + 1)
meta_log <- log2(assay(metabric_common_h5, i = assay_meta)[common_genes, ] + 1)

var_tcga <- rowVars(tcga_log)
var_meta <- rowVars(meta_log)

mean_var <- (var_tcga + var_meta) / 2
names(mean_var) <- common_genes

candidate_genes <- names(sort(mean_var, decreasing = TRUE))[1:2000]


# 16. Objective 3: Gene-wise Cox 

run_univariate_cox <- function(se, assay_name = NULL, genes, surv_df) {
  common_ids <- intersect(colnames(se), surv_df$sample_id)
  surv_df <- surv_df[match(common_ids, surv_df$sample_id), ]
  
  out <- data.frame(
    gene = genes,
    beta = NA_real_,
    HR = NA_real_,
    p = NA_real_
  )
  
  for (i in seq_along(genes)) {
    g <- genes[i]
    
    expr <- if (is.null(assay_name)) {
      as.numeric(log2(assay(se)[g, common_ids] + 1))
    } else {
      as.numeric(log2(assay(se, i = assay_name)[g, common_ids] + 1))
    }
    
    if (sd(expr, na.rm = TRUE) == 0) next
    
    expr_z <- as.numeric(scale(expr))
    
    fit <- try(
      coxph(Surv(time_years, event) ~ expr_z, data = surv_df),
      silent = TRUE
    )
    
    if (!inherits(fit, "try-error")) {
      s <- summary(fit)
      out$beta[i] <- s$coefficients[1, "coef"]
      out$HR[i]   <- s$coefficients[1, "exp(coef)"]
      out$p[i]    <- s$coefficients[1, "Pr(>|z|)"]
    }
  }
  
  out$q <- p.adjust(out$p, method = "BH")
  out
}

res_tcga <- run_univariate_cox(
  se = tcga_common_h5,
  assay_name = assay_tcga,
  genes = candidate_genes,
  surv_df = surv_tcga
)

res_meta <- run_univariate_cox(
  se = metabric_common_h5,
  assay_name = assay_meta,
  genes = candidate_genes,
  surv_df = surv_meta
)

# 17. Objective 3: Filter and Identify Shared Genes

res_tcga <- res_tcga[
  is.finite(res_tcga$beta) &
    is.finite(res_tcga$HR) &
    !is.na(res_tcga$p) &
    res_tcga$HR > 0 &
    res_tcga$HR < 10,
]

res_meta <- res_meta[
  is.finite(res_meta$beta) &
    is.finite(res_meta$HR) &
    !is.na(res_meta$p) &
    res_meta$HR > 0 &
    res_meta$HR < 10,
]

sig_tcga <- res_tcga[
  !is.na(res_tcga$q) & res_tcga$q < 0.05,
]

sig_meta <- res_meta[
  !is.na(res_meta$q) & res_meta$q < 0.05,
]

merged_sig <- merge(
  sig_tcga,
  sig_meta,
  by = "gene",
  suffixes = c("_tcga", "_meta")
)

shared_candidate <- merged_sig[
  merged_sig$beta_tcga * merged_sig$beta_meta > 0,
]

shared_candidate$rank_score <-
  -log10(shared_candidate$q_tcga) +
  -log10(shared_candidate$q_meta)

shared_candidate <- shared_candidate[
  order(shared_candidate$rank_score, decreasing = TRUE),
]


# 18. Objective 3: Build Candidate Signature 

candidate_signature <- shared_candidate$gene

cand_score_tcga <- score_signature(tcga_common_h5, assay_tcga, candidate_signature)
cand_score_meta <- score_signature(metabric_common_h5, assay_meta, candidate_signature)

surv_tcga$candidate_score <- cand_score_tcga[surv_tcga$sample_id]
surv_meta$candidate_score <- cand_score_meta[surv_meta$sample_id]

cox_tcga <- coxph(Surv(time_years, event) ~ candidate_score, data = surv_tcga)
cox_meta <- coxph(Surv(time_years, event) ~ candidate_score, data = surv_meta)


# 19. Objective 3: KM Plots

surv_meta <- surv_meta[!is.na(surv_meta$candidate_score), ]
surv_meta$sig_group <- ifelse(
  surv_meta$candidate_score >= median(surv_meta$candidate_score, na.rm = TRUE),
  "High signature", "Low signature"
)

fit_meta <- survfit(Surv(time_years, event) ~ sig_group, data = surv_meta)

ggsurvplot(
  fit_meta,
  data = surv_meta,
  palette = c("#0072B2", "#CC79A7"),
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (years)",
  ylab = "Overall survival probability",
  title = "METABRIC survival stratified by candidate gene signature",
  legend.title = "Candidate signature",
  legend.labs = c("High", "Low")
)

surv_tcga <- surv_tcga[!is.na(surv_tcga$candidate_score), ]
surv_tcga$sig_group <- ifelse(
  surv_tcga$candidate_score >= median(surv_tcga$candidate_score, na.rm = TRUE),
  "High signature", "Low signature"
)

fit_tcga <- survfit(Surv(time_years, event) ~ sig_group, data = surv_tcga)

ggsurvplot(
  fit_tcga,
  data = surv_tcga,
  palette = c("#0072B2", "#CC79A7"),
  risk.table = TRUE,
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (years)",
  ylab = "Overall survival probability",
  title = "TCGA-BRCA survival stratified by candidate gene signature",
  legend.title = "Candidate signature",
  legend.labs = c("High", "Low")
)


# 20. Objective 3: HR Comparison Plot

plot_data <- shared_candidate[, c("gene", "HR_tcga", "HR_meta")]
plot_data$gene <- factor(plot_data$gene, levels = rev(plot_data$gene))

ggplot(plot_data, aes(x = gene)) +
  geom_point(aes(y = HR_tcga, shape = "TCGA"), size = 3) +
  geom_point(aes(y = HR_meta, shape = "METABRIC"), size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_classic(base_size = 12) +
  labs(
    title = "Cross-cohort prognostic effects of candidate genes",
    x = "Gene",
    y = "Hazard ratio",
    shape = "Cohort"
  )


# 21. Objective 3: Tables 

# Shared candidate genes table
table_genes <- shared_candidate[, c(
  "gene",
  "HR_tcga", "q_tcga",
  "HR_meta", "q_meta",
  "rank_score"
)]

table_genes[, -1] <- round(table_genes[, -1], 3)
write.csv(table_genes, "objective3_shared_candidate_genes.csv", row.names = FALSE)

# Signature performance table
table_signature <- rbind(
  extract_cox_summary(cox_tcga, "TCGA-BRCA"),
  extract_cox_summary(cox_meta, "METABRIC")
)

table_signature[, -1] <- round(table_signature[, -1], 3)
write.csv(table_signature, "objective3_signature_performance.csv", row.names = FALSE)


# 22.Dataset summary table


# Total expression samples
tcga_total_samples <- ncol(tcga_common_h5)
meta_total_samples <- ncol(metabric_common_h5)

# Shared genes after harmonisation
shared_genes_n <- length(intersect(rownames(tcga_common_h5), rownames(metabric_common_h5)))

# Samples used in PAM50 survival analysis
tcga_surv_samples <- nrow(surv_tcga_pam50)
meta_surv_samples <- nrow(surv_meta_pam50)

dataset_summary <- data.frame(
  Dataset = c("TCGA-BRCA", "METABRIC"),
  `Data type` = c("RNA sequencing", "Microarray"),
  Platform = c("Illumina HiSeq", "Illumina HT-12"),
  `Expression data` = c("Raw counts (log2-transformed)", "Microarray intensity (log2-transformed)"),
  `Total expression samples` = c(tcga_total_samples, meta_total_samples),
  `Samples included in PAM50 survival analysis` = c(tcga_surv_samples, meta_surv_samples),
  `Shared genes after harmonisation` = c(shared_genes_n, shared_genes_n),
  check.names = FALSE
)

dataset_summary
write.csv(dataset_summary, "dataset_summary_table.csv", row.names = FALSE)


