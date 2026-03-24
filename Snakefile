import os
from pathlib import Path
import pandas as pd

rebuild_tags = config.get("rebuild_tags", False)

df = pd.read_csv('data/cohort.csv', dtype=str).fillna("")

def norm(x: str) -> str:
    return (x or "").strip()

# ----------------------------
# Metadata
# ----------------------------
meta = {}
samples_all = []
samples_vis = []
samples_scrna = []
samples_with_probes = []

samples_exemplar = set([
    'Visium_HD_Human_Breast_Cancer_FFPE',
    'Visium_HD_Human_Kidney_FFPE',
    'Visium_HD_Human_Lung_Cancer_Fixed_Frozen',
    'Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1',
    'Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2',
    'Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1',
    'Visium_HD_Human_Lymph_Node_FFPE',
    'Visium_HD_Human_Prostate_Cancer_FFPE',
    'Visium_HD_Human_Tonsil_Fresh_Frozen',
    'Visium_HD_Human_Tonsil_Fresh_Frozen_IF', 'Visium_HD_Mouse_Brain',
    'Visium_HD_Mouse_Brain_Fixed_Frozen',
    'CytAssist_11mm_FFPE_Human_Kidney',
    'CytAssist_11mm_FFPE_Human_Glioblastoma',
    'SC3_v3_NextGem_DI_PBMC_CSP_1K', 'Parent_NGSC3_DI_PBMC',
    '5k_pbmc_v3_nextgem', 'manual_5k_pbmc_NGSC3_ch5',
    'pbmc_10k_protein_v3', '1k_PBMCs_TotalSeq_B_3p', 'pbmc_10k_v3',
    'SC3_v3_NextGem_DI_PBMC_10K', 'pbmc_1k_v3',
    'connect_5k_pbmc_NGSC3_ch1', '5k_pbmc_v3',
    'connect_5k_pbmc_NGSC3_ch5', 'pbmc_1k_protein_v3',
    'SC3_v3_NextGem_SI_PBMC_CSP_1K', 'manual_5k_pbmc_NGSC3_ch1',
    'SC3_v3_NextGem_SI_PBMC_10K', '10k_PBMC_3p_nextgem_Chromium_X',
    'pbmc1k_v2', 'pbmc4k_v2', 'pbmc8k_v2',
    '20k_PBMC_5pv2_HT_nextgem_Chromium_X', 'Lung_Cancer_PBMC',
    'PBMC_Control2', 'Kidney_Sarcomatoid_PBMC', 'Kidney_Cancer_PBMC',
    'PBMC_Control1',
    '10k_Human_PBMC_TotalSeqB_singleplex_10k_Human_PBMC_TotalSeqB_singleplex',
    '16plex_human_PBMC_TotalSeqC_multiplex_PBMC_01_BC1_AB1',
    'PBMC_F1B_10X_3p_v31', 'PBMC_F5B_10X_3p_v31',
    'PBMC_F1A_10X_3p_v31', 'PBMC_F5A_10X_3p_v31',
    'PBMC_F5A_10X_flex_v1', 'PBMC_F1B_10X_flex_v1',
    'PBMC_F1A_10X_flex_v1', 'PBMC_F5B_10X_flex_v1'])

samples_for_overlaps = set([
    'Visium_HD_Human_Breast_Cancer_FFPE','Visium_HD_Human_Kidney_FFPE',
    'Visium_HD_Human_Lung_Cancer_HD_Only_Experiment1',
    'Visium_HD_Human_Lung_Cancer_post_Xenium_Prime_5K_Experiment2',
    'Visium_HD_Human_Lung_Cancer_post_Xenium_v1_Experiment1',
    'Visium_HD_Human_Lymph_Node_FFPE', 'Visium_HD_Human_Prostate_Cancer_FFPE',
    'Visium_HD_Human_Tonsil_Fresh_Frozen',
    'Visium_HD_Human_Tonsil_Fresh_Frozen_IF', 'Visium_HD_Mouse_Brain',
    'PBMC_F1B_10X_3p_v31', 'PBMC_F5B_10X_3p_v31', 'PBMC_F1A_10X_3p_v31',
    'PBMC_F5A_10X_3p_v31','PBMC_F5A_10X_flex_v1', 'PBMC_F1B_10X_flex_v1',
    'PBMC_F1A_10X_flex_v1','PBMC_F5B_10X_flex_v1',
    'SC3_v3_NextGem_DI_PBMC_CSP_1K', 'Parent_NGSC3_DI_PBMC',
    '5k_pbmc_v3_nextgem', 'manual_5k_pbmc_NGSC3_ch5',
    'pbmc_10k_protein_v3', 'pbmc_10k_v3', 'SC3_v3_NextGem_DI_PBMC_10K',
    'pbmc_1k_v3', 'connect_5k_pbmc_NGSC3_ch1', '5k_pbmc_v3',
    'pbmc_1k_protein_v3', 'SC3_v3_NextGem_SI_PBMC_CSP_1K',
    'manual_5k_pbmc_NGSC3_ch1', 'SC3_v3_NextGem_SI_PBMC_10K',
    'Lung_Cancer_PBMC', 'PBMC_Control2', 'Kidney_Sarcomatoid_PBMC',
    'Kidney_Cancer_PBMC', 'PBMC_Control1',
    '10k_Human_PBMC_TotalSeqB_singleplex_10k_Human_PBMC_TotalSeqB_singleplex'])


for _, r in df.iterrows():
    dataset = norm(r["dataset"])
    sample = norm(r["sample"])
    if not dataset or not sample:
        continue

    if not rebuild_tags and not os.path.isfile(os.path.join(dataset, sample, 'scdepth_tags.gz')):
        print('Skipping ', dataset, sample)


    bam = norm(r["bam"])
    gtf = norm(r["gtf"])
    lib = norm(r["library_prep"])
    positions = norm(r.get("positions", ""))
    probes = norm(r.get("probes", ""))

    # ext inferred from positions filepath (only used for visium/visium_hd targets)
    pos_ext = ""
    if positions:
        suf = Path(positions).suffix.lower().lstrip(".")
        pos_ext = suf if suf in {"csv", "parquet"} else "csv"

    meta[(dataset, sample)] = {
        "bam": bam,
        "gtf": gtf,
        "lib": lib,
        "positions": positions,
        "pos_ext": pos_ext,
        "probes": probes,
    }

    samples_all.append((dataset, sample))

    if dataset in {"visium", "visium_hd"}:
        samples_vis.append((dataset, sample))

    if dataset in {"scrna", "scrna_flex", "scrna_expanded"}:
        samples_scrna.append((dataset, sample))

    if probes:
        samples_with_probes.append((dataset, sample))

tdf = pd.read_csv('data/visium_tissue_data.csv', keep_default_na=False)
tmap = {k.sample:(k.tissue_frac, k.tissue_type) for k in tdf.itertuples() if k.tissue_frac > 0.0}
CACHE_TAGS = []
if rebuild_tags:
    CACHE_TAGS = [f"{d}/{s}/scdepth_tags.gz" for d, s in samples_all]
CACHE_BINDEX = [f"{d}/{s}/scdepth_barcode_index.txt.gz" for d, s in samples_all]
CACHE_SUMMARY = [f"{d}/{s}/scdepth_summary.txt" for d, s in samples_all]

VISIUM_POS_OUT = []
for d, s in samples_vis:
    if meta[(d, s)]["positions"]:
        ext = meta[(d, s)]["pos_ext"] or "csv"
        VISIUM_POS_OUT.append(f"{d}/{s}/scdepth_positions.{ext}")

EXCLUDE_OUT = [f"{d}/{s}/scdepth_exclude.txt" for d, s in samples_with_probes]

EMPTYDROPS_OUT = [f"{d}/{s}/scdepth_emptydrops.txt.gz" for d, s in samples_scrna]

FIT_BASELINE = [f"{d}/{s}/scdepth_fit_baseline.txt" for d, s in samples_all]
FIT_CURVE = [f"{d}/{s}/scdepth_fit_curve.txt" for d, s in samples_all]
FIT_SVG = [f"{d}/{s}/scdepth_fit_fits.svg" for d, s in samples_all]

LIMIT_SUMMARY = [f"{d}/{s}/scdepth_limit_summary.txt" for d, s in samples_all if s in samples_exemplar]
PRESEQ_SUMMARY = [f"{d}/{s}/scdepth_preseq_summary.txt" for d, s in samples_all if s in samples_exemplar]
GO_SUMMARY = [f"{d}/{s}/scdepth_genes_summary.txt" for d, s in samples_all if s in samples_for_overlaps]
ST_SUMMARY = [f"{d}/{s}/scdepth_stability_summary.txt" for d, s in samples_all if s in samples_for_overlaps]

rule all:
    input:
        # cache outputs
        CACHE_TAGS,
        CACHE_BINDEX,
        CACHE_SUMMARY,

        # optional per-dataset outputs
        VISIUM_POS_OUT,
        EXCLUDE_OUT,
        EMPTYDROPS_OUT,

        # fits
        FIT_BASELINE,
        FIT_CURVE,
        FIT_SVG,

        #limit analyses
        LIMIT_SUMMARY,
        PRESEQ_SUMMARY,

	#gene overalaps
	GO_SUMMARY,
	ST_SUMMARY,


rule cache_tags:
    input:
        gtf=lambda wc: meta[(wc.dataset, wc.sample)]["gtf"],
        bam=lambda wc: meta[(wc.dataset, wc.sample)]["bam"],
    output:
        tags=protected("{dataset}/{sample}/scdepth_tags.gz"),
        bindex=protected("{dataset}/{sample}/scdepth_barcode_index.txt.gz"),
        summary=protected("{dataset}/{sample}/scdepth_summary.txt"),
    params:
        lib=lambda wc: meta[(wc.dataset, wc.sample)]["lib"],
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 10
    run:
        # If tags already exists, skip scdepth cache entirely and do not modify anything.
        if os.path.exists(output.tags):
            missing = [p for p in [output.bindex, output.summary] if not os.path.exists(p)]
            if missing:
                raise RuntimeError(
                    f"{output.tags} exists but other cache outputs are missing: {missing}. "
                    "Refusing to re-run cache to avoid touching existing tags."
                )
            print(f"[cache_tags] Found existing {output.tags}; skipping cache.")
            return

        shell("scdepth cache -t 4 -l {params.lib} -g {input.gtf} {input.bam} {params.prefix}")

rule copy_barcodes:
    input:
        lambda wc: meta[(wc.dataset, wc.sample)]["positions"]
    output:
        "{dataset}/{sample}/scdepth_positions.{ext}"
    wildcard_constraints:
        dataset="visium|visium_hd",
        ext="csv|parquet"
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 30
    shell:
        "scdepth barcodes -b {input} {params.prefix}"

rule exclude_probes:
    input:
        lambda wc: meta[(wc.dataset, wc.sample)]["probes"]
    output:
        "{dataset}/{sample}/scdepth_exclude.txt"
    wildcard_constraints:
        dataset="visium|visium_hd|scrna_flex"
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 30
    shell:
        "scdepth probes -p {input} {params.prefix}"

rule scrna_emptydrops:
    input:
        tags="{dataset}/{sample}/scdepth_tags.gz",
        bindex="{dataset}/{sample}/scdepth_barcode_index.txt.gz",
        summary="{dataset}/{sample}/scdepth_summary.txt",
    output:
        cells="{dataset}/{sample}/scdepth_emptydrops.txt.gz",
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    wildcard_constraints:
        dataset="scrna|scrna_flex|scrna_expanded"
    priority: 20
    shell:
        "scdepth emptydrops -t 4 {params.prefix}"

def fit_inputs(wc):
    base = [
        f"{wc.dataset}/{wc.sample}/scdepth_tags.gz",
        f"{wc.dataset}/{wc.sample}/scdepth_barcode_index.txt.gz",
        f"{wc.dataset}/{wc.sample}/scdepth_summary.txt",
    ]
    if meta[(wc.dataset, wc.sample)]["probes"]:
        base.append(f"{wc.dataset}/{wc.sample}/scdepth_exclude.txt")
    return base

rule build_curve:
    input:
        fit_inputs
    output:
        "{dataset}/{sample}/scdepth_fit_baseline.txt",
        "{dataset}/{sample}/scdepth_fit_curve.txt",
        "{dataset}/{sample}/scdepth_fit_fits.svg",
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 100
    shell:
        "scdepth fit -s {wildcards.sample} -t 5 {params.prefix}"

def pilot_data(wc):
    if 'scrna' in wc.dataset:
        return '--use-scrna'
    if wc.sample in tmap:
        t = tmap[wc.sample]
        if t[1] != '':
            return f'--tissue-frac {t[0]} --capture-format {t[1]}'
        else:
            return f'--tissue-frac {t[0]}'

    return ''

def limits_inputs(wc):
    ins = fit_inputs(wc)[:]  # tags/bindex/summary (+ exclude if probes)
    if wc.dataset in {"scrna", "scrna_flex", "scrna_expanded"}:
        ins.append(f"{wc.dataset}/{wc.sample}/scdepth_emptydrops.txt.gz")
    return ins

def pos_inputs(wc):
    ins = limits_inputs(wc)
    if wc.dataset in {'visium', 'visium_hd'}:
        ext = meta[(wc.dataset, wc.sample)]["pos_ext"] or 'csv'
        ins.append(f"{wc.dataset}/{wc.sample}/scdepth_positions.{ext}")
    return ins


rule run_limits:
    input:
        #fit_inputs,
        limits_inputs,
        baseline="{dataset}/{sample}/scdepth_fit_baseline.txt",
    output:
        "{dataset}/{sample}/scdepth_limit_summary.txt",
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
        pilot=pilot_data
    priority: 150
    shell:
        "scdepth limit -s {wildcards.sample} {params.pilot} -t 5 {params.prefix}"

rule run_preseq:
    input:
        "{dataset}/{sample}/scdepth_limit_summary.txt",
    output:
        "{dataset}/{sample}/scdepth_preseq_summary.txt",
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 100
    shell:
        "scdepth preseq -s {wildcards.sample} --preseq=preseq/build/preseq {params.prefix}"

rule run_goverlaps:
    input:
        pos_inputs,
        baseline="{dataset}/{sample}/scdepth_fit_baseline.txt",
    output:
        "{dataset}/{sample}/scdepth_genes_summary.txt",
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 100
    shell:
        "scdepth genes -s {wildcards.sample} -t 5 {params.prefix}"

rule run_stab:
    input:
        pos_inputs,
        baseline="{dataset}/{sample}/scdepth_fit_baseline.txt",
    output:
        "{dataset}/{sample}/scdepth_stability_summary.txt",
    resources:
        parallel_limit=3
    params:
        prefix=lambda wc: f"{wc.dataset}/{wc.sample}/scdepth",
    priority: 100
    shell:
        "scdepth stability -s {wildcards.sample} -t 5 {params.prefix}"
