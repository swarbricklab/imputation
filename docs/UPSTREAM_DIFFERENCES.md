# Differences from the upstream sc-eQTLGen pipeline

This workflow is derived from the **sc-eQTLGen consortium WG1 genotype QC + imputation
pipeline** (Powell Lab / Garvan; maintainer Drew Neavin), the same group behind Demuxafy.

- Upstream repo: <https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC>
- Upstream docs: <https://wg1-pipeline-qc.readthedocs.io/>
- Key scripts we inherited: `Imputation/scripts/filter_het.R`,
  `PCA_Projection_Plotting.R` (author header: *Jose Alquicira Hernandez, 2019*).

Upstream ships as a single monolithic Singularity image
(`SNP_imputation_1000g_hg38.sif`) built from manually-assembled binaries, and was
designed for **germline PBMC/blood** cohorts (e.g. OneK1K). This document records how
our version differs, so the changes are auditable for publication and so the deltas are
explicit when comparing outputs.

This pipeline reproduces the upstream/`.sif` ("golden") output **byte-for-byte except for
cosmetic header timestamps and record/sample ordering** — see *Reproduction status* below.

The differences below were verified rule-by-rule against the vanilla upstream
(`sc-eQTLgen-consortium/WG1-pipeline-QC@a9a81bd`,
`Imputation/includes/{urmo_imputation_hg38.smk,plink_gender_ancestry_QC.smk}`). Every
divergence is deliberate; §2.1–2.6 are changes made here, while §2.8 records customisations
already baked into the retired `.sif` (they reproduce the golden but differ from vanilla
upstream).

---

## 1. Technical differences (packaging & reproducibility)

| Area | Upstream | This version |
|------|----------|--------------|
| **Containerisation** | one monolithic `.sif` of hand-assembled binaries | per-rule `conda:` environments (`workflow/envs/*.yaml`) plus a few pinned biocontainer images for `bcftools`. Tool versions are pinned to match the `.sif` (CrossMap 0.6.5, plink 1.90b6.21 / plink2 2.00a3.7, GenotypeHarmonizer 1.4.23, bcftools 1.10.2, vcftools 0.1.16, Eagle 2.4.1, Minimac4 1.0.2). |
| **Reference / tool data** | baked into the `.sif`, or fetched by hand | tracked publicly via `dvc import-url` and unpacked by `prep.sh` (the eQTLGen imputation reference bundle, 1000G, the GRCh37→GRCh38 chain, GenotypeHarmonizer, the Minimac4 static binary). Transparent, versioned provenance. |
| **Configurability** | several thresholds hard-coded in the rules | exposed as config params: `pre_maf`, `post_maf`, `snp_missing_pct`, `snp_hwe`, `chunk_length`, `het_remove_outliers`. |
| **Determinism** | final VCF carries run-varying headers (Minimac `##filedate`, bcftools `; Date=` stamps and `##bcftools_*Command` lines with per-run temp paths); intermediate `temp.vcf` left on disk | `rename_chromosomes` strips all run-varying header metadata, so reruns from identical inputs are **byte-identical**; `temp.vcf` is marked `temp()`. |
| **Chromosome ordering** | final VCF is string-sorted (`chr1, chr10, chr11, …, chr2, …`) | naturally ordered (`chr1..chr22, chrX`) via `bcftools sort` over a chr-named `--fai` reheader. |
| **Snakemake compatibility** | older Snakemake | runs under Snakemake 7.32.4; the `het_filter` param was renamed `remove` → `remove_outliers` because `remove` is reserved and otherwise fails to parse. |
| **Execution** | site-specific | portable Snakemake profiles: `profiles/global/nci` (PBS Pro fan-out) and `profiles/workflow` (per-rule resources); runs as a git submodule inside a DVC super-project. |

---

## 2. Logical / processing differences

### 2.1 Heterozygosity-rate sample filter (`scripts/filter_het.R`) — **bug fix**
Upstream classifies samples by heterozygosity rate `(N_SITES - O(HOM)) / N_SITES` and
removes those outside ±3 SD. The "passed" subset is built with a **malformed boolean**:

```r
# upstream (buggy): the two OR clauses collapse to just (HET_RATE < mean + 3*sd)
het_pass <- subset(het, (HET_RATE < mean-3*sd) | (HET_RATE < mean+3*sd))
```

Because the condition reduces to `HET_RATE < mean + 3*sd`, upstream only ever removes
**high**-heterozygosity outliers and **silently keeps every low-heterozygosity outlier**
(a sample can even appear in *both* the `failed` and `passed` lists). Upstream then
removes the complement unconditionally (`bcftools view -S passed_list`).

Our fix uses the symmetric band and makes removal opt-in:

```r
# fixed: symmetric ±3 SD
het_pass <- subset(het, (HET_RATE > mean-3*sd) & (HET_RATE < mean+3*sd))
```

- `het_filter` gains a `het_remove_outliers` param (**default `false`** → flag outliers in
  the report but keep all samples). Set `true` to restore removing behaviour.
- **Rationale:** tumour-derived genotypes show low heterozygosity from loss of
  heterozygosity (LOH); dropping them discards good donors and is harmful for
  demultiplexing. Low het is expected signal here, not contamination.
- **Reproduction-neutral:** with `het_remove_outliers: false` and no high-het outliers
  present, this keeps exactly the same sample set the upstream bug did (so it does not
  change the golden output for this cohort). The upstream bug has not yet been reported.

### 2.2 X/Y reinsertion (`get_preimpute_XY`, `reinsert_XY`) — **robustness fix**
Autosomes are imputed per-ancestry; X/Y are carried through from the pre-imputation VCF
and merged back. Upstream assumes the autosome and X/Y sample sets are identical, which
crashes whenever upstream QC (or the het filter) drops any sample. Our version:
- `get_preimpute_XY` subsets X/Y to the post-QC sample set;
- `reinsert_XY` derives the sample order from the imputed autosomes and subsets X/Y to
  match before `bcftools concat`.

This fixes the crash and makes the merge order well-defined. (Note: the merged sample
*column order* therefore differs from upstream — a data-preserving permutation.)

### 2.3 Ancestry PCA assignment (`PCA_Projection_Plotting.R`) — **bug fix**
Upstream clusters with `pam(scores[, 2:11], …)`, which positionally includes the IID
column instead of ten principal components; under R ≥ 4 this collapses borderline
ancestry calls. Our version selects the PC columns explicitly (`pam(scores[, pc_cols], …)`).

**Automated sex/ancestry decisions.** Upstream writes sex-mismatch and ancestry-mismatch
rows with a blank `UPDATE/REMOVE/KEEP` column and **halts for a human to curate them
interactively**. Our version is fully automated: it sets `UPDATE` for every sex mismatch
(sex is corrected to the genotype-inferred value) and every ancestry mismatch (ancestry is
set to the PCA assignment). Combined with the auto-generated psam (§2.7; `SEX=0`,
`Provided_Ancestry=NONE`), this means every sample's sex and ancestry are assigned from the
genotype/PCA, and **no sample is dropped for a sex or ancestry mismatch** (only the
relatedness check, §2.4, removes samples). This is a deliberate de-interactivation for
unattended HPC runs — appropriate for an array cohort with no curated phenotype, but note it
is not a provided-vs-inferred QC gate.

### 2.4 Relatedness checking — **new feature**
A `relatedness_check` rule was added to the QC stage to flag related sample pairs
(absent upstream).

### 2.5 Sample-metadata / ancestry outputs — **new feature**
The QC stage exports per-sample metadata (updated sex / ancestry assignments,
`export_sample_metadata`) as first-class outputs.

### 2.6 Post-imputation MAF default
Upstream (and the retired `.sif`) hard-coded `MAF >= 0.05` in the post-imputation
`filter_maf_r2` step. We expose it as `post_maf` and **currently default it to `0.05`**
to reproduce the golden output. A more permissive `0.01` is planned but is held back
until a publishable version that reproduces the golden exists.

### 2.7 Array/VCF input front-end — **new**
Upstream starts from a fully-populated plink2 pgen/pvar/psam (real `SEX` and
`Provided_Ancestry`). This version accepts the genotyping stage's VCF and adapts it:
- `reheader_vcf` relabels sample IDs from `id_map`;
- `run_plink` converts VCF → pgen (`--make-pgen --sort-vars`), attaching a psam generated by
  `scripts/generate_sample_psam.py` with `SEX=0` and `Provided_Ancestry='NONE'` for all
  samples (the array supplies no phenotype);
- `restore_vcf_header` (imputation stage) renames samples back to the original donor IDs via
  `id_map` before the final output.

Consequence: because the generated psam carries no provided sex/ancestry, the QC checks
always register a mismatch, which §2.3's auto-`UPDATE` then fills from the genotype/PCA. This
is intended for an array cohort with no supplied phenotype.

### 2.8 Post-imputation filtering — customisations **inherited from the retired `.sif`**
These were already baked into the `.sif` / golden pipeline (not introduced here): they
reproduce the golden output but differ from the vanilla sc-eQTLGen
`urmo_imputation_hg38.smk`. Documented because they materially change the variant set:
- **`filter_maf_r2` info-score threshold `R2 > 0.8`** — vanilla upstream uses `R2 >= 0.3`; a
  much stricter imputation-quality cut.
- **`filter_maf_r2` filters imputed variants only** —
  `(IMPUTED=1 && MAF >= {post_maf} && R2 > 0.8) || (IMPUTED=0)` keeps all genotyped/typed
  (`IMPUTED=0`) variants unconditionally; vanilla upstream applies the MAF/R² cut to every
  variant.
- **No `bcftools +fill-tags`** immediately before `filter_maf_r2` — the filter uses the MAF
  tag Minimac4 emits rather than recomputing MAF across the cohort.
- **No `--max-missing 1` "complete_cases" step** — vanilla upstream's `filter4demultiplexing`
  drops every variant with any missing genotype after the exon/indel filter; ours retains
  them.
- **`filter_preimpute_vcf` uses a single global `pre_maf`** for all ancestries, vs vanilla
  upstream's per-ancestry MAF table. Moot at the configured `pre_maf: 0.0` (no pre-imputation
  MAF filtering), but a real change in where the threshold comes from.

Both `R2 > 0.8` and the imputed-only filter are **intended** (confirmed): they are the
behaviour the golden `.sif` produced, which this pipeline reproduces byte-for-byte.

---

## 3. Reproduction status

Running this version (conda/biocontainers, `post_maf: 0.05`, `het_remove_outliers: false`)
on identical inputs reproduces the golden `.sif` output **byte-for-byte** after
normalising cosmetic differences:

- identical sample set (211) and variant set (215,875);
- every autosome and chrX byte-identical per-chromosome (records **and** INFO; genotype
  and dosage differences zero);
- whole-genome records md5 identical once sample columns and records are sorted.

The only residual differences are **(a)** header timestamps (now stripped for
determinism), **(b)** natural vs string-sorted chromosome order, and **(c)** sample
column order (from §2.2). Eagle 2.4.1 and Minimac4 1.0.2 produced bit-identical
phasing/imputation across fan-out nodes — there is no numerical drift.
