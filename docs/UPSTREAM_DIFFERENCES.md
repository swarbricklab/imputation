# Differences from the upstream sc-eQTLGen WG1 pipeline

This workflow is adapted from the **sc-eQTLGen consortium WG1 genotype-QC + imputation
pipeline** (Powell Lab / Garvan; maintainer Drew Neavin — the group behind
[Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/)). The WG1 pipeline
was built for germline PBMC/blood cohorts to support the OneK1K project; our adaptations are
partly intended to support tumour cohorts, where some loss of heterozygosity (LOH) is to be
expected.

- Upstream repo: <https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC>
  (pinned reference for the comparison below: `@a9a81bd`,
  `Imputation/includes/{urmo_imputation_hg38.smk, plink_gender_ancestry_QC.smk}`)
- Upstream docs: <https://wg1-pipeline-qc.readthedocs.io/>
- Inherited scripts: `scripts/filter_het.R`, `PCA_Projection_Plotting.R`
  (*Jose Alquicira Hernandez, 2019*).

> **A note on lineage.** The Powell Lab also publishes
> [`powellgenomicslab/SNP_imputation_1000g_hg38`](https://github.com/powellgenomicslab/SNP_imputation_1000g_hg38),
> a **sibling** of WG1 — *not* our parent and not WG1's ancestor. The two repos share a
> common sc-eQTLGen origin and diverged on 2021-11-16 (when the Powell repo was spun out);
> Powell then switched ancestry assignment from `pam` clustering to `knn` (2022-06) and
> carries no post-imputation filtering. This workflow is on the **WG1 line** — it keeps
> `pam` clustering and an R²/MAF/exon post-imputation filter, both of which Powell dropped
> after the split. The table includes Powell to make the distinction explicit; everything
> after the table compares against WG1.

## At a glance

| Aspect | **This workflow** | **WG1-pipeline-QC** `@a9a81bd` | Powell `SNP_imputation_1000g_hg38` `@437b0bd` |
|---|---|---|---|
| Input | VCF → pgen (`reheader_vcf`; generated psam `SEX=0`/`NONE`) | starts from pgen (real `SEX`/`Provided_Ancestry`) | starts from pgen |
| het `het_pass` logic | symmetric `&` (fixed) | **OR-bug** — keeps low-het outliers | symmetric `&` |
| het outlier removal | opt-out, **default keep all** | always removes | always removes |
| Ancestry method | `pam(PC1–10)`, explicit cols | `pam(scores[,2:11])` — **IID-col bug** | `knn` (caret, PC1–10) |
| Sex/ancestry mismatch | auto-`UPDATE`, unattended | blank col, **interactive** curation | blank col, interactive |
| Pre-imputation MAF | single global `pre_maf` (config) | per-ancestry, **interactive prompt** | per-ancestry, interactive |
| Chromosomes | 1–22 **+ X/Y reinserted** | 1–22 only (X/Y dropped) | 1–22 only |
| Post-imputation filter | `filter_maf_r2` + `filter_exons_indels` | `filter4demultiplexing` | **none** |
| &nbsp;&nbsp;· R² | **`> 0.8`**, configurable | `>= 0.3`, hardcoded | — |
| &nbsp;&nbsp;· MAF | `>= post_maf`, configurable (0.05) | `>= 0.05`, hardcoded | — |
| &nbsp;&nbsp;· scope | imputed-only (keeps typed) | all variants | — |
| &nbsp;&nbsp;· `+fill-tags` first | no (uses Minimac MAF) | yes | — |
| &nbsp;&nbsp;· complete-cases (`--max-missing 1`) | **no** | yes | — |
| Final chr naming | `rename_chromosomes`: natural order, headers stripped (deterministic) | string order, run-varying | string order, run-varying |
| Relatedness check | yes | no | no |
| Sample-metadata export | yes | no | no |
| Tooling | per-rule public OCI images (+conda fallback) | monolithic `.sif` | monolithic `.sif` |
| Reference data | public `dvc import-url` + `prep.sh` | baked / manual | baked / manual |

The sections below explain the *why* behind these differences (all vs WG1).

## Processing differences

### Bug fixes

**Heterozygosity-rate filter** (`scripts/filter_het.R`). WG1's "passed" subset uses a
malformed boolean that collapses to `HET_RATE < mean + 3*sd` (both clauses use `<`), so it
only removes *high*-het outliers and silently keeps every *low*-het one:

```r
# WG1 (buggy): both OR clauses reduce to (HET_RATE < mean + 3*sd)
het_pass <- subset(het, (HET_RATE < mean-3*sd) | (HET_RATE < mean+3*sd))
# ours: symmetric +/-3 SD
het_pass <- subset(het, (HET_RATE > mean-3*sd) & (HET_RATE < mean+3*sd))
```

We use the symmetric band and make removal opt-in via `het_remove_outliers` (**default
`false`**: flag outliers in the report, keep all samples). Rationale: tumour-derived
genotypes have low heterozygosity from LOH; dropping them discards good donors and harms
demultiplexing. (The sibling Powell repo independently uses the symmetric band too.)

**Ancestry-PCA clustering** (`PCA_Projection_Plotting.R`). WG1 clusters with
`pam(scores[, 2:11], …)`, which positionally includes the IID column instead of ten PCs,
collapsing borderline ancestry calls under R ≥ 4. We select the PC columns explicitly
(`pam(scores[, pc_cols], …)`).

### Behavioural changes

**Fully automated, unattended operation.** WG1 prompts the operator interactively — for a
per-ancestry pre-imputation MAF, and to curate sex/ancestry mismatches (it writes a blank
`UPDATE/REMOVE/KEEP` column and halts). This version runs unattended: it auto-sets `UPDATE`
for every sex and ancestry mismatch (to the genotype/PCA value) and takes thresholds from
config, so **no sample is dropped for a sex/ancestry mismatch** (only the relatedness check
removes samples). Appropriate for an array cohort with no curated phenotype; it is not a
provided-vs-inferred QC gate.

**X/Y reinsertion** (`get_preimpute_XY`, `reinsert_XY`). WG1 images autosomes 1–22 only and
discards X/Y. We carry X/Y through from the pre-imputation VCF, subset them to the post-QC
sample set, and merge them back (order derived from the imputed autosomes). Consequence:
merged sample *column order* differs from the autosome-only output — a data-preserving
permutation.

**Array/VCF input front-end.** WG1 starts from a populated plink2 pgen/pvar/psam (real `SEX`,
`Provided_Ancestry`). This version takes the genotyping stage's VCF: `reheader_vcf` relabels
sample IDs from `id_map`; `run_plink` builds a pgen with a generated psam (`SEX=0`,
`Provided_Ancestry=NONE`); `restore_vcf_header` renames samples back before the final output.
With no provided phenotype, every sample registers a QC mismatch that the auto-`UPDATE` fills.

**Single global pre-imputation MAF.** WG1 applies a per-ancestry MAF (from the interactive
prompt); we use one global `pre_maf` config value (default `0.0` — no pre-imputation MAF
filter). Both run `bcftools +fill-tags` first.

### New outputs

- **Relatedness check** (`relatedness_check`) — flags related sample pairs (absent in WG1).
- **Sample-metadata export** (`export_sample_metadata`) — per-sample updated sex/ancestry as
  first-class outputs.

### Post-imputation filtering

WG1's `filter4demultiplexing` runs `+fill-tags` → `bcftools filter 'MAF>=0.05 & R2>=0.3'`
(every variant) → exon/indel filter → `--max-missing 1` complete-cases. Our `filter_maf_r2`
+ `filter_exons_indels` keep the same *shape* — filter by MAF and R², then restrict to exons
— but differ, trading variant count for confidence (fewer, higher-quality SNPs for
demultiplexing):

- **R² `> 0.8`** (WG1: `>= 0.3`) and **MAF `>= post_maf`** — both now **config params**
  (`post_r2`, `post_maf`; defaults `0.8` / `0.05`), so a site can retune without editing rules.
- **imputed variants only**: `(IMPUTED=1 && MAF>={post_maf} && R2>0.8) || (IMPUTED=0)` keeps
  all genotyped (typed) variants unconditionally; WG1 applies its MAF/R² cut to every variant.
- **no `+fill-tags`** before the filter — uses Minimac4's emitted MAF rather than recomputing
  it across the cohort.
- **no `--max-missing 1`** complete-cases step — variants with any missing genotype are kept.
- **`rename_chromosomes`** then imposes natural `chr1..chr22,chrX` order and strips
  run-varying header metadata, so reruns from identical inputs are byte-identical (WG1 leaves
  string-sorted output with run-varying headers). The sibling Powell repo has no
  post-imputation filtering at all.

## Packaging & reproducibility

- **Containers, not a monolith.** Each rule runs a pinned public tool — a per-rule OCI image
  by default, or conda as a fallback — replacing the retired monolithic Singularity image.
  Tool versions are pinned to those the retired `.sif` shipped (bcftools 1.10.2, plink
  1.90b6.21 / plink2 2.00a3.7, vcftools 0.1.16, CrossMap 0.6.5, GenotypeHarmonizer 1.4.23,
  Eagle 2.4.1, Minimac4 1.0.2). See the README.
- **Public reference data** via `dvc import-url` + `prep.sh`, vs baked-in / manual refs.
- **Configurable thresholds** (`pre_maf`, `post_maf`, `post_r2`, `snp_missing_pct`, `snp_hwe`,
  `chunk_length`, `het_remove_outliers`) vs hard-coded / interactive.
- **Deterministic & unattended** — no interactive prompts; run-varying header metadata is
  stripped, so reruns from identical inputs are byte-identical; `temp.vcf` is `temp()`.
- **Portable execution** — Snakemake 7.32.4 (the `het_filter` param `remove` → `remove_outliers`,
  since `remove` is reserved); `nci` (PBS) + workflow profiles; runs as a submodule in a DVC
  super-project.
