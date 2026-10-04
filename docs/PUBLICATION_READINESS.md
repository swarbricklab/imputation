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

> **Status (2026-10).** Most of Goal 1 is now done: references are public
> (#22/#27), tool versions are recorded and every rule runs from a public
> per-rule OCI image (#20/#23 via #31/#33 — done ahead of schedule), the profile
> and preflight are harmonised with genotyping (#26 via #32), and the rule graphs
> are regenerated against the current rules. **The remaining Goal-1 work is
> citability (§1): make the repo public, review `CITATION.cff`'s author list, then
> tag v1.0.0 + wire Zenodo for the DOI.** Two small loose ends remain in §4
> (`dvc push` the bundle; retire the now-unused `.sif` import). Items below are
> annotated inline.

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
- [x] Regenerate the rule graphs and confirm they match the code (issue #11).
      Regenerated 2026-10 against the current rules (QC 18 rules, imputation 22);
      the previous SVGs were from 2024 and had drifted (18 nodes). NB: like
      genotyping, the super-project's `run_*.sh` writes graphs to the
      super-project `docs/`, not this repo's `docs/`, so they can drift —
      regenerate deliberately after rule changes.

### 3. Record & pin tool versions (issue #20)
- [x] Tool versions recorded in `docs/UPSTREAM_DIFFERENCES.md` §1 (CrossMap 0.6.5,
      plink 1.90b6.21 / plink2 2.00a3.7, GenotypeHarmonizer 1.4.23, bcftools
      1.10.2, vcftools 0.1.16, Eagle 2.4.1, Minimac4 1.0.2), pinned to match the
      retired `.sif`. Full graduation to per-tool **public** images was completed
      ahead of schedule (#31/#33 — see §6), so this exceeds the Goal-1 minimum.

### 4. Public reference data (issue #22) — the main Goal-1 reproducibility blocker
The QC and imputation references are private `dvc import`s from
`Swarbricklab/references.git`. But those are a **repackaging of a public
bundle**: the sceQTL-Gen consortium / Powell Lab distribute the whole imputation
reference set as `eQTLGenImpRef.tar.gz` (see the upstream
[wiki](https://github.com/powellgenomicslab/SNP_imputation_1000g_hg38/wiki/SNP-Genotype-Imputation-Using-1000G-hg38-Reference)).
So most of the private refs can be re-sourced from **one** public download.

**Verified-live public sources** (checked 2026-09):

- **`eQTLGenImpRef.tar.gz`** — `https://www.dropbox.com/s/l60a2r3e4vo78mn/eQTLGenImpRef.tar.gz?dl=1`
  (md5: `https://www.dropbox.com/s/eci808v0uepqgcz/eQTLGenImpRef.tar.gz.md5?dl=1`).
  Unpacks to `hg38/` with `imputation/` (Minimac4 ref), `phasing/genetic_map/`
  (Eagle maps), `phasing/phasing_reference/`, `ref_genome_QC/` (GRCh38 primary
  assembly), `ref_panel_QC/` (the 30x-GRCh38 panel). This one bundle supplies
  **five** of the refs below.
- The monolithic container `SNP_imputation_1000g_hg38.sif` is also public —
  `https://www.dropbox.com/s/mwjpndpclp6njcg/SNP_imputation_1000g_hg38.sif?dl=1`
  (relevant to #20/#23: pin/record it from this source).

| Config key | Current (private) | Public source |
|---|---|---|
| `refs.vcf` (imputation ref panel) | `…/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted.vcf.gz` | **eQTLGenImpRef bundle** → `hg38/ref_panel_QC/` |
| `refs.hg38_int_fa` | `…/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa` | **eQTLGenImpRef bundle** → `hg38/ref_genome_QC/` |
| `refs.genetic_map` | `…/hg38/phasing/genetic_map/…` | **eQTLGenImpRef bundle** → `hg38/phasing/genetic_map/` |
| `refs.phasing` | `…/hg38/phasing/phasing_reference/` | **eQTLGenImpRef bundle** → `hg38/phasing/phasing_reference/` |
| `refs.impute` | `resources/reference/hg38/imputation` | **eQTLGenImpRef bundle** → `hg38/imputation/` |
| `refs.hg38_chr_fai` | `…/refdata-gex-GRCh38-2020-A/fasta/genome.fa.fai` | prefer the bundle's `ref_genome_QC` FASTA `.fai` (or genotyping's public GRCh38) — **verify contigs match what the panel expects** |
| `refs.1000g.*` (QC ancestry) | `resources/1000g/all_phase3_filtered.*` | 1000G phase-3 plink (public; document the exact filtering) — **not in the bundle** |
| `refs.bed` | `resources/bed/hg38exonsUCSC.bed` | UCSC exon BED (public `dvc import-url`) |
| `refs.*chr_map*` | `resources/genomes/chr_map/*` | small; ship in-repo or regenerate |

**Recipe** (belongs in the super-project, where the data lives — see the brca rewire):
```bash
# fetch on a data-mover node (internet-capable queue)
qx exec --internet --env dt3 -P a56 --storage gdata/a56+scratch/a56 \
  -- dvc import-url 'https://www.dropbox.com/s/l60a2r3e4vo78mn/eQTLGenImpRef.tar.gz?dl=1' \
     resources/imputation/eQTLGenImpRef.tar.gz
# then a dvc stage extracts hg38/ into the paths the config expects
```
- [x] Bundle **tracked via `dvc import-url`** from the public source (fetched on
      a data-mover node with `qx --internet`, ~36 GB): `resources/eQTLGenImpRef.tar.gz.dvc`.
      Out md5 `88e3603933a21712a7403023c1f2e9df` matches the published checksum —
      byte-identical to the canonical public reference, so equivalence to the old
      private import holds by construction (the private import was extracted from
      this same bundle).
- [x] **`prep.sh`** extracts the bundle into the config paths
      (`resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}`,
      `resources/reference/hg38/imputation`).
- [x] Retired the two private imports the bundle covers: `resources/reference.dvc`
      (79 GB Minimac4 ref) and `resources/genomes/hg38.dvc` (44 GB panel/phasing/QC).
- [ ] `dvc push` the bundle to the imputation remote (data-mover). *(outstanding)*
- [x] **`bed` and `1000g` are now public** — both tracked via public `dvc
      import-url` (Dropbox); `chr_map` ships in-repo. All reference imports are
      public *except* the vestigial `resources/imputation/SNP_imputation_1000g_hg38.sif.dvc`.
- [ ] **Retire `resources/imputation/SNP_imputation_1000g_hg38.sif.dvc`** — it is
      the only remaining *private* import and is no longer used by any rule (the
      workflow graduated off the monolithic `.sif` in #33; it survives only in
      comments). Remove it so the reference set is fully public. *(outstanding)*

### 5. Harmonisation with genotyping (shared front-end; issue #26 in genotyping)
- [x] Switch the nested `profiles/global` submodule URL **ssh → https** (public).
- [x] `run_qc.sh` / `run_imputation.sh` now use `profiles/global/nci` (the
      generalised, env-templated `PROJECT` / `jobmode=pbspro` profile), not
      `nci_a56`; submodule pinned accordingly.
- [x] Adopted genotyping's `prep.sh` read-only **preflight** (`--what check`):
      implemented in `prep.sh` and called from `run_qc.sh` / `run_imputation.sh`
      with a `SKIP_PREFLIGHT` opt-out.

---

## Goal 2 — improvements (after v1.0.0)

### 6. Graduate from the monolithic container (issue #23) — ✅ DONE (ahead of schedule)
Completed in #31/#33: `SNP_imputation_1000g_hg38.sif` is fully retired; every rule
now runs from a per-rule **public** OCI image on `ghcr.io/swarbricklab`
(`bcftools`, `plink`/`plink-bcftools`, `vcftools`, `minimac4`, `eagle`, `crossmap`,
`het`, and a hand-built **GenotypeHarmonizer** image — no public biocontainer
existed). Images are referenced via a `config['containers']` map (centralised in
#35) with pinned tags; `conda:` directives are retained but commented so a
conda-only fallback stays one edit away. The last cleanup is retiring the unused
`.sif` DVC import (see §4).

### 7. Keep imputation permissive; move the tunable threshold to snp_demux (#24)
`rule vcf_filter` (postimputation.smk) hardcodes
`bcftools filter --include 'MAF>=0.05 & R2>=0.3'` **and** bakes the thresholds
into output filenames. Change to:
- carry `INFO/R2` (and MAF) through unfiltered, or apply only a permissive floor;
- parameterise `R2`/`MAF` in config (no thresholds in filenames);
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
