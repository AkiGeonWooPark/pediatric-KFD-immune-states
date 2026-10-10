#!/usr/bin/env python3
"""
Single-cell TCR repertoire analysis (scirpy v0.22.1)

`adata_gex` is the annotated KFD PBMC AnnData (exported from Seurat);
filtered_contig_annotations.csv files come from Cell Ranger multi (V(D)J-T).
"""
import numpy as np
import anndata as ad
import muon as mu
import scirpy as ir

# ---- 1. Load contigs and pair with matched GEX ------------------------------
tcr = ad.concat([ir.io.read_10x_vdj(f) for f in contig_files],
                join="outer", merge="same")
mdata = mu.MuData({"gex": adata_gex, "airr": tcr})
common = mdata["gex"].obs_names.intersection(mdata["airr"].obs_names)
mdata = mdata[common, :]

# ---- 2. Chain indexing + QC -------------------------------------------------
ir.pp.index_chains(mdata)
ir.tl.chain_qc(mdata)
mu.pp.filter_obs(mdata, "airr:chain_pairing", lambda x: x != "multichain")
mu.pp.filter_obs(mdata, "airr:chain_pairing",
                 lambda x: ~np.isin(x, ["orphan VDJ", "orphan VJ"]))
# restrict to conventional alpha-beta T cells
mdata = mdata[(mdata.obs["gex:Broad_Anno"] == "T")
              & ~mdata.obs["gex:Anno1"].isin(["rd_T", "NK"])].copy()

# ---- 3. TCRdist clonotype clusters (same V gene) -----------------------------
ir.pp.ir_dist(mdata, metric="tcrdist", sequence="aa", cutoff=15)
ir.tl.define_clonotype_clusters(
    mdata, sequence="aa", metric="tcrdist",
    receptor_arms="all", dual_ir="any",
    same_v_gene=True, key_added="cc_aa_tcrdist_same_v",
)

# ---- 4. Clonal expansion ----------------------------------------------------
# Expanded = clonotype cluster containing >= 2 cells.
ir.tl.clonal_expansion(mdata, target_col="cc_aa_tcrdist_same_v")
mdata.obs[["airr:cc_aa_tcrdist_same_v", "airr:cc_aa_tcrdist_same_v_size",
           "airr:clonal_expansion"]].to_csv("TCR_clonotypes.csv")
