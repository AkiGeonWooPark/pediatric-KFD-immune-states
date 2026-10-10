# =============================================================================
# Single-cell RNA-seq preprocessing

# PBMC gene-expression, hashtag (TotalSeq-C) and V(D)J libraries were processed
# with `cellranger multi` (v7.1.0) against GRCh38-2020-A (GEX) and
# GRCh38-alts-ensembl-7.1.0 (V(D)J). Three KFD captures (GEM1-3) and one
# healthy-control capture were processed identically; this script shows the
# per-capture flow, followed by per-sample QC.
# =============================================================================

library(Seurat)          # v5
library(SoupX)
library(scDblFinder)
library(SingleCellExperiment)

# ---- 1. Ambient RNA correction (SoupX, per capture) -------------------------
run_soupx <- function(raw, filt, force_rho = NULL) {
  tmp <- CreateSeuratObject(filt$`Gene Expression`, min.cells = 3, min.features = 200)
  tmp <- NormalizeData(tmp) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA()
  tmp <- FindNeighbors(tmp, dims = 1:10) |> FindClusters(resolution = 0.5) |> RunUMAP(dims = 1:10)
  genes <- intersect(rownames(raw$`Gene Expression`), rownames(filt$`Gene Expression`))
  sc <- SoupChannel(raw$`Gene Expression`[genes, ],
                    filt$`Gene Expression`[genes, colnames(tmp)])
  sc <- setClusters(sc, setNames(tmp$seurat_clusters, colnames(tmp)))
  sc <- autoEstCont(sc, forceAccept = TRUE)
  if (!is.null(force_rho)) sc$metaData$rho <- force_rho
  adjustCounts(sc)
}
counts <- run_soupx(raw, filt)

# ---- 2. Hashtag demultiplexing ----------------------------------------------
obj <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
obj[["HTO"]] <- CreateAssayObject(filt$`Antibody Capture`[, colnames(obj)])
obj <- NormalizeData(obj, assay = "HTO", normalization.method = "CLR")
obj <- HTODemux(obj, assay = "HTO", positive.quantile = 0.99)
obj <- subset(obj, subset = HTO_classification.global == "Singlet")

# ---- 3. Per-sample adaptive QC ----------------------------------------------
qc_filter <- function(s, nmads = 3) {
  s[["percent.mt"]] <- PercentageFeatureSet(s, pattern = "^MT-")
  s$log10GenesPerUMI <- log10(s$nFeature_RNA) / log10(s$nCount_RNA)
  m <- function(x) median(x); d <- function(x) mad(x)
  min_f  <- max(300, m(s$nFeature_RNA) - nmads * d(s$nFeature_RNA))
  max_f  <-          m(s$nFeature_RNA) + nmads * d(s$nFeature_RNA)
  min_c  <- max(500, m(s$nCount_RNA)   - nmads * d(s$nCount_RNA))
  max_mt <- min(15,  m(s$percent.mt)   + nmads * d(s$percent.mt))
  min_cx <- max(0.80, m(s$log10GenesPerUMI) - nmads * d(s$log10GenesPerUMI))
  subset(s, subset = nFeature_RNA > min_f & nFeature_RNA < max_f &
                     nCount_RNA > min_c & percent.mt < max_mt &
                     log10GenesPerUMI > min_cx)
}
samples <- lapply(samples, qc_filter)

# ---- 4. Doublet removal (scDblFinder, per sample) ---------------------------
samples <- lapply(samples, function(s) {
  sce <- scDblFinder(as.SingleCellExperiment(s))
  s$scDblFinder.class <- sce$scDblFinder.class
  subset(s, subset = scDblFinder.class == "singlet")
})

# ---- 5. Merge ---------------------------------------------------------------
# Two donors were captured in two runs (GEM1 and GEM2); their samples are
# merged per donor, then all donors are merged into one object.
obj <- merge(samples[[1]], y = samples[-1])
saveRDS(obj, "pbmc_qc.rds")
