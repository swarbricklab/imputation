# Publication-readiness plan

This document tracks the work to make `imputation` publication-ready, following
the same two-goal framing used for the
[genotyping](https://github.com/swarbricklab/genotyping) workflow:

- **Goal 1 — publication-ready *without changing processing*.** Make the repo
  public, citable, and reproducible from public inputs, without altering the
  numbers the workflow produces. This must be finished first.
- **Goal 2 — improvements.** Changes that alter behaviour or structure
  (container graduation, permissive R², parameterisation, etc.).

The genotyping workflow completed Goal 1 and released v1.0.0; this plan reuses
the patterns established there.

---

## Goal 1 — publication-ready (no processing change)

### 1. Citability (issue #19)
- [x] `LICENSE` (MIT, Garvan Institute of Medical Research) — mirrors genotyping.
- [x] `CITATION.cff` — **draft**; the author list is copied from genotyping and
      **must be reviewed** for imputation-specific contributors before release.
- [ ] Make the repository **public**.
- [ ] Enable the Zenodo–GitHub webhook, tag **v1.0.0**, add the concept DOI to
      `CITATION.cff` (`identifiers:`). Same route as genotyping (webhook, not the
      manual reserve path).

### 2. Documentation & provenance (issue #21)
- [x] README rewrite: fixed typos/broken links, documented the two workflows,
      inputs/outputs, tools, and reference-data provenance + citations.
- [ ] Confirm the citation lineage with the authors (Võsa/eQTLGen →
      sceQTL-Gen WG1 → powellgenomicslab/SNP_imputation_1000g_hg38 (Neavin) →
      swarbricklab; van der Wijst 2020 eLife; Demuxafy / Neavin 2024 Genome Biol).
- [ ] Regenerate the rule graphs once the workflow is runnable from public
      inputs, and confirm they match the code (issue #11). NB: like genotyping,
      the super-project's `run_*.sh` writes graphs to the super-project `docs/`,
      not this repo's `docs/`, so they can drift — regenerate deliberately.

### 3. Record & pin tool versions (issue #20)
- [ ] The monolithic `SNP_imputation_1000g_hg38.sif` is opaque. As a Goal-1
      *minimum*, record the exact version of every tool it bundles (bcftools,
      plink/plink2, vcftools, Minimac4, Eagle, GenotypeHarmonizer, CrossMap,
      Java, R + tidyverse/ggpubr/cluster/RColorBrewer/dplyr) and its build
      provenance, so runs are reproducible and documented.
      (Full graduation to per-tool images is Goal 2 / issue #23.)

### 4. Public reference data (issue #22) — the main Goal-1 reproducibility blocker
The QC and imputation references are private `dvc import`s from
`Swarbricklab/references.git` (no public URL). This is the exact problem
genotyping solved for GRCh38 (rebuild from a public source via `dvc import-url`,
verify byte-/variant-equivalence). Replace each with a public source:

| Config key | Current (private) | Public source to use |
|---|---|---|
| `refs.1000g.*` (QC ancestry) | `resources/1000g/all_phase3_filtered.*` | 1000G phase 3 plink files (public; document the filtering) |
| `refs.vcf` (imputation ref panel) | `resources/genomes/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted.vcf.gz` | 1000G **30x GRCh38** panel (public) |
| `refs.hg38_int_fa` | `.../ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | Ensembl GRCh38 primary assembly (as genotyping did) |
| `refs.hg38_chr_fai` | `resources/genomes/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai` | **Same private 10x genome genotyping just retired** — reuse genotyping's public GRCh38 build |
| `refs.genetic_map` | `.../phasing/genetic_map/genetic_map_hg38_withX.txt.gz` | Eagle genetic maps (public) |
| `refs.phasing` | `.../phasing/phasing_reference/` | 1000G phasing reference (public) |
| `refs.impute` | `resources/reference/hg38/imputation` | Minimac4 1000G reference (public) |
| `refs.bed` | `resources/bed/hg38exonsUCSC.bed` | UCSC exon BED (public) |
| `refs.*chr_map*` | `resources/genomes/chr_map/*` | small, ship in-repo or regenerate |
- [ ] Verify variant-level equivalence after each swap (the technique used for
      genotyping: reconstruct/compare, ignore header-only differences).

### 5. Harmonisation with genotyping (shared front-end; issue #26 in genotyping)
- [x] Switch the nested `profiles/global` submodule URL **ssh → https** (public).
- [ ] Update the `snakemake_config` pin to the **generalised `nci` profile**
      (env-templated `PROJECT`, `jobmode=pbspro`) that genotyping v1.0.0 uses,
      and update `run_qc.sh` / `run_imputation.sh` from `profiles/global/nci_a56`
      to `profiles/global/nci`; require `PROJECT` in the environment.
- [ ] Adopt genotyping's `prep.sh` + read-only **preflight** (`--what check`)
      pattern to verify containers/references/inputs before submitting jobs, and
      call it from the run scripts (with a `SKIP_PREFLIGHT` opt-out).

---

## Goal 2 — improvements (after v1.0.0)

### 6. Graduate from the monolithic container (issue #23)
Replace `SNP_imputation_1000g_hg38.sif` (used by ~29 rule invocations) with
per-rule **public** biocontainers: `bcftools` (already partly done), `plink2`,
`vcftools`, `minimac4`, `eagle`, `crossmap`, and a **GenotypeHarmonizer** image
published to GHCR/Docker Hub (no public biocontainer exists — build one, as
discussed). Pin every tag.

### 7. Keep imputation permissive; move the tunable threshold to snp_demux (#24)
`rule filter_maf_r2` (imputation.smk) hardcodes `R2 > 0.8` on imputed sites
(MAF is already set by `post_maf`). Change to:
- carry `INFO/R2` (and MAF) through unfiltered, or apply only a permissive floor;
- parameterise `R2` in config;
- move the operational, tunable R² cut to `snp_demux` (see snp_demux #37), so a
  demux parameter sweep does not require re-running imputation.

### 8. Other
- [ ] Extract ancestry mappings (issue #15).
- [ ] Installation script (issue #5).
- [ ] Revisit the QC/`plink_qc` boundary with genotyping (which intrinsic QC, if
      any, belongs upstream — cf. genotyping #31). Decision so far: **keep
      ancestry here**; evaluate only the SNP-level QC steps.

---

## Notes carried over from the genotyping work
- Use the Zenodo **webhook** route (repo must be public first; can't pre-reserve
  into that snapshot).
- Large public downloads run on a data-mover node (`qx --internet --env dt3`).
- Reference/data reuse: prefer copying content-addressed objects between DVC
  caches over re-downloading (as done for the mega-atlas GRCh38).
- The `dt` wrapper isn't on the default PATH on all nodes; git hooks may fail
  with `dt: not found` — commits/pushes there need `--no-verify`.
