# =============================================================================
# Type I IFN module score and severity-associated cell-type signatures
# =============================================================================

library(Seurat)
library(msigdbr)
library(dplyr)

# ---- 1. Hallmark IFN-alpha response module score ----------------------------
ifna <- msigdbr(species = "Homo sapiens", category = "H") |>
  filter(gs_id == "M5911") |> pull(gene_symbol) |> unique()
obj <- AddModuleScore(obj, features = list(ifna), name = "H_IFNA_Res")

# ---- 2. Signature gene sets for spatial projection --------------------------
Idents(obj) <- "Anno1"
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sig <- markers |>
  filter(p_val_adj < 0.05,
         (cluster %in% c("pDC", "CD8_T_Memory_3") & avg_log2FC > 2) |
         (cluster %in% c("Classical_Mono_IFN_stim",
                         "Non_Classical_Mono_IFN_stim") & avg_log2FC > 4))
signatures <- split(sig$gene, droplevels(sig$cluster))

# Specificity check in the scRNA-seq data
obj <- AddModuleScore(obj, features = signatures, name = "signature")
saveRDS(signatures, "signature_gene_sets.rds")
