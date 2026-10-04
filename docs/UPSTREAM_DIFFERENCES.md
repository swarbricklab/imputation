# Differences from the upstream sc-eQTLGen pipeline

This workflow is adapted from the **sc-eQTLGen consortium WG1 genotype-QC + imputation
pipeline** (Powell Lab / Garvan; maintainer Drew Neavin — the group behind Demuxafy),
which ships as one monolithic Singularity image (`SNP_imputation_1000g_hg38.sif`) built
for germline PBMC/blood cohorts (e.g. OneK1K).

- Upstream repo: <https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC>
- Upstream docs: <https://wg1-pipeline-qc.readthedocs.io/>
- Inherited scripts: `scripts/filter_het.R`, `PCA_Projection_Plotting.R`
  (*Jose Alquicira Hernandez, 2019*).

Every divergence below is deliberate and was verified rule-by-rule against vanilla
upstream (`sc-eQTLgen-consortium/WG1-pipeline-QC@a9a81bd`).

## Does it change the output? No.

On identical inputs this pipeline reproduces the retired `.sif` ("golden") output
**byte-for-byte** after normalising cosmetic differences:

- identical sample set (211) and variant set (215,875);
- every autosome and chrX byte-identical per chromosome — records, INFO, genotypes and
  dosages, zero differences;
- whole-genome records md5-identical once sample columns and records are sorted.

The only residual differences are cosmetic: **(a)** header timestamps (deliberately
stripped for determinism), **(b)** natural `chr1..chr22,chrX` vs upstream's string-sorted
chromosome order, and **(c)** sample-column order (from the X/Y fix below). Eagle 2.4.1 and
Minimac4 1.0.2 produce bit-identical phasing/imputation across nodes — no numerical drift.

## Processing differences

### Bug fixes

**Heterozygosity-rate filter** (`scripts/filter_het.R`). Upstream's "passed" subset uses a
malformed boolean that collapses to `HET_RATE < mean + 3*sd`, so it only removes *high*-het
outliers and silently keeps every *low*-het one:

```r
# upstream (buggy): both OR clauses reduce to (HET_RATE < mean + 3*sd)
het_pass <- subset(het, (HET_RATE < mean-3*sd) | (HET_RATE < mean+3*sd))
# fixed: symmetric +/-3 SD
het_pass <- subset(het, (HET_RATE > mean-3*sd) & (HET_RATE < mean+3*sd))
```

We fix the band and make removal opt-in via `het_remove_outliers` (**default `false`**: flag
outliers in the report, keep all samples). Rationale: tumour-derived genotypes have low
heterozygosity from loss of heterozygosity (LOH); dropping them discards good donors and
harms demultiplexing. Reproduction-neutral — with the default and no high-het outliers, the
retained set is exactly what the upstream bug produced.

**Ancestry-PCA clustering** (`PCA_Projection_Plotting.R`). Upstream's `pam(scores[, 2:11], …)`
positionally includes the IID column instead of ten PCs, which collapses borderline ancestry
calls under R ≥ 4. We select the PC columns explicitly.

### Behavioural changes

**Fully automated sex/ancestry assignment.** Upstream writes mismatch rows with a blank
`UPDATE/REMOVE/KEEP` column and halts for interactive curation. This version auto-sets
`UPDATE` for every sex and ancestry mismatch (correcting to the genotype/PCA value), so it
runs unattended and **drops no sample for a sex/ancestry mismatch** (only the relatedness
check removes samples). Appropriate for an array cohort with no curated phenotype; it is not
a provided-vs-inferred QC gate.

**X/Y reinsertion** (`get_preimpute_XY`, `reinsert_XY`). Autosomes are imputed per ancestry
and X/Y carried through from the pre-imputation VCF. Upstream assumes identical autosome and
X/Y sample sets and crashes whenever QC drops a sample; we subset X/Y to the post-QC set and
derive the merge order from the imputed autosomes. Consequence: merged sample *column order*
differs from upstream — a data-preserving permutation.

**Array/VCF input front-end.** Upstream starts from a populated plink2 pgen/pvar/psam (real
`SEX`, `Provided_Ancestry`). This version takes the genotyping stage's VCF: `reheader_vcf`
relabels sample IDs from `id_map`; `run_plink` builds a pgen with a generated psam
(`SEX=0`, `Provided_Ancestry=NONE`); `restore_vcf_header` renames samples back before the
final output. Because the psam carries no provided phenotype, every sample registers a QC
mismatch that the auto-`UPDATE` above fills from the genotype/PCA.

### New outputs

- **Relatedness check** (`relatedness_check`) — flags related sample pairs (absent upstream).
- **Sample-metadata export** (`export_sample_metadata`) — per-sample updated sex/ancestry as
  first-class outputs.

### Post-imputation filtering

The `filter_maf_r2` step carries customisations **baked into the golden `.sif`** that differ
from vanilla upstream (`urmo_imputation_hg38.smk`) and materially change the variant set:

- **`R2 > 0.8`** info-score cut (vanilla: `R2 >= 0.3`) — much stricter.
- **imputed variants only**: `(IMPUTED=1 && MAF >= {post_maf} && R2 > 0.8) || (IMPUTED=0)`
  keeps all typed variants unconditionally (vanilla filters every variant).
- **no `bcftools +fill-tags`** beforehand — uses Minimac4's emitted MAF, not a recomputed one.
- **no `--max-missing 1`** complete-cases step — variants with missing genotypes are retained
  (vanilla's `filter4demultiplexing` drops them).
- **single global `pre_maf`** for all ancestries (vanilla: a per-ancestry MAF table); moot at
  the configured `pre_maf: 0.0`.

`post_maf` is exposed (default **`0.05`** to reproduce the golden; a permissive `0.01` is
planned post-release). `R2 > 0.8` and the imputed-only filter are confirmed intended.

## Packaging & reproducibility

- **Containers, not a monolith.** Each rule runs a pinned public tool — a per-rule OCI image
  by default, or conda as a fallback — replacing the opaque `.sif`. Versions match the `.sif`
  (bcftools 1.10.2, plink 1.90b6.21 / plink2 2.00a3.7, vcftools 0.1.16, CrossMap 0.6.5,
  GenotypeHarmonizer 1.4.23, Eagle 2.4.1, Minimac4 1.0.2). See the README.
- **Public reference data** via `dvc import-url` + `prep.sh`, vs hand-assembled/baked refs.
- **Configurable thresholds** (`pre_maf`, `post_maf`, `post_r2`, `snp_missing_pct`, `snp_hwe`,
  `chunk_length`, `het_remove_outliers`) vs hard-coded.
- **Deterministic output** — run-varying header metadata (Minimac `##filedate`, bcftools
  `Date=` / command lines) is stripped, so reruns are byte-identical; `temp.vcf` is `temp()`.
- **Natural chromosome order** (`chr1..chr22,chrX`) vs upstream's string sort.
- **Snakemake 7.32.4** (the `het_filter` param `remove` → `remove_outliers`, since `remove`
  is reserved); portable `nci` (PBS) + workflow profiles; runs as a submodule in a DVC
  super-project.
