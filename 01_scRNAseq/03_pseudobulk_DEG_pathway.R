# =============================================================================
# Pseudobulk differential expression (KFD vs HC) and pathway enrichment
# =============================================================================

library(Seurat)
library(clusterProfiler)
library(org.Hs.eg.db)

# ---- 1. Pseudobulk --------------------
obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
pb <- AggregateExpression(obj, assays = "RNA",
                          group.by = c("disease", "GEM", "donor", "Broad_Anno"),
                          return.seurat = TRUE)
pb$group <- paste(pb$Broad_Anno, pb$disease, sep = "_")
Idents(pb) <- "group"

# ---- 2. Batch-aware gene filter + Wilcoxon DE (per cell type) ---------------
de_celltype <- function(ct) {
  cols  <- pb$Broad_Anno == ct
  mat   <- GetAssayData(pb, layer = "counts")[, cols]
  batch <- pb$GEM[cols]
  keep  <- Reduce(`&`, lapply(unique(batch), function(b)
             rowSums(mat[, batch == b, drop = FALSE] > 0) > 0))
  de <- FindMarkers(pb, ident.1 = paste0(ct, "_KFD"), ident.2 = paste0(ct, "_HC"),
                    features = rownames(mat)[keep])          # Wilcoxon (default)
  de$gene <- rownames(de)
  de
}
de <- lapply(setNames(nm = c("T", "B", "Mo_DC")), de_celltype)

# DEGs: P < 0.05 and |log2FC| > 1

# ---- 3. GO / KEGG enrichment ------------------------------------------------
enrich <- function(genes) {
  ids <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
  list(GO   = enrichGO(ids, OrgDb = org.Hs.eg.db, ont = "BP",
                       pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE),
       KEGG = enrichKEGG(ids, organism = "hsa", pvalueCutoff = 0.05))
}
