#!/usr/bin/env python3
"""
Cell-cell interaction inference (CellPhoneDB v5.0.1, database v5.0.0)

Inputs exported from the annotated Seurat object:
  counts_norm.txt  - library-size-normalized counts (counts / total * 1e4)
  meta.tsv         - barcode -> annotated cell type
"""
from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

cpdb_results = cpdb_statistical_analysis_method.call(
    cpdb_file_path="cellphonedb.zip",       # database v5.0.0
    meta_file_path="meta.tsv",
    counts_file_path="counts_norm.txt",
    counts_data="hgnc_symbol",
    iterations=1000,                        # label permutations
    threshold=0.1,                          # min. fraction of expressing cells
    pvalue=0.05,                            # permutation P (no multiplicity correction)
    score_interactions=True,
    debug_seed=42,
    threads=5,
    output_path="cellphonedb_out",
)
