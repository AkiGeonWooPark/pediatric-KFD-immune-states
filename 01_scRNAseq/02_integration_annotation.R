# =============================================================================
# Integration, clustering and annotation (Seurat v5 + Harmony)
# =============================================================================

library(Seurat)          # v5
library(harmony)

run_flow <- function(obj, resolution) {
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj, selection.method = "vst")
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  obj <- RunHarmony(obj, group.by.vars = "donor")      # default parameters
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
  FindClusters(obj, resolution = resolution)           # Louvain (default)
}

# ---- 1. Broad clustering ----------------------------------------------------
obj <- run_flow(obj, resolution = 0.2)
markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# ---- 2. Lineage re-clustering -----------------------------------------------
t_nk    <- run_flow(subset(obj, idents = t_nk_clusters),    resolution = 1.0)
b       <- run_flow(subset(obj, idents = b_clusters),       resolution = 1.0)
myeloid <- run_flow(subset(obj, idents = myeloid_clusters), resolution = 1.0)

# ---- 3. Merge annotated lineages --------------------------------------------
obj <- merge(pdc, y = list(t_nk, b, myeloid))
obj <- NormalizeData(obj) |> FindVariableFeatures() |> ScaleData() |> RunPCA()
obj <- RunHarmony(obj, group.by.vars = "donor")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
saveRDS(obj, "pbmc_annotated.rds")

