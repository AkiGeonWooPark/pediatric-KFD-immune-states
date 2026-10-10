# =============================================================================
# GeoMx DSP (Human Whole Transcriptome Atlas) analysis
#
# Segment/probe QC, target filtering, background correction and Q3
# normalization were performed in the GeoMx DSP Data Analysis Suite
# (v3.1.0); the exported Q3-normalized target-count matrix is the input here.
# =============================================================================

library(dplyr)
library(purrr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(escape)

# ---- 1. ROI-level differential expression -----------------------------------
deg <- map_dfr(genes, function(g) {
  x <- log2(expr[[g]] + 1)
  tibble(Gene = g,
         logFC = mean(x[expr$severity_group == "high"]) -
                 mean(x[expr$severity_group == "low"]),
         p_value = wilcox.test(x ~ expr$severity_group)$p.value)
}) |> mutate(FDR = p.adjust(p_value, method = "BH"))
up <- deg |> filter(FDR < 0.05, logFC > 0.5) |> pull(Gene)


# ---- 2. Pathway enrichment of severity-high genes ---------------------------
ids  <- bitr(up, "SYMBOL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
go   <- enrichGO(ids, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")
kegg <- enrichKEGG(ids, organism = "hsa")

# ---- 3. Projection of PBMC-derived signatures (AUCell via escape) -----------
signatures <- readRDS("signature_gene_sets.rds")   # from 01_scRNAseq/04
auc <- escape.matrix(t(as.matrix(expr[, genes])), gene.sets = signatures,
                     method = "AUCell", groups = 30, min.size = NULL)
# AUCell scores compared between severity groups (Wilcoxon) and correlated
# with severity score (Spearman).
