# imputation

A Snakemake workflow that takes **measured genotypes** (a VCF, e.g. from the
[genotyping](https://github.com/swarbricklab/genotyping) workflow) and prepares
them for single-cell demultiplexing by quality-controlling the calls and
**imputing** missing variants against a 1000 Genomes reference panel.

It contains two workflows that run in sequence:

1. **QC** (`Snakefile_qc`, `rules/qc.smk`) — filters SNPs, and checks reported
   sex and assigns ancestry by comparison against a 1000 Genomes reference.
2. **Imputation** (`Snakefile_imputation`, `rules/imputation.smk`,
   `rules/postimputation.smk`) — lifts the genotypes to GRCh38, harmonises and
   fixes reference alleles, phases, imputes missing variants against a 1000
   Genomes reference panel, and filters the result (MAF, imputation R², exonic
   regions) to produce SNP profiles.

These are the second and third stages of the SNP-based demultiplexing pipeline:

```
genotyping  ->  imputation (this repo: QC + imputation)  ->  snp_demux (Demuxafy)
```

The final VCF of imputed SNP profiles is consumed by
[Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/) to
assign cells to donors in multiplexed 10x Chromium pools.

## Inputs (from `genotyping`)

- **VCF** — genotype calls for all samples (hg19 for QC, hg38 for imputation)
- **psam** — a sample description file for `plink`
- **id map** — maps original sample ids to the `plink`-safe ids in the `.psam`

## Outputs

- **QC**: a post-QC set of `plink` files (`pgen`/`pvar`/`psam`) with sex
  double-checked and ancestry assigned from the 1000 Genomes reference.
- **Imputation**: imputed, filtered SNP profiles ready for Demuxafy.

## Rule graphs

![QC rule graph](docs/qc_rulegraph.svg)

![Imputation rule graph](docs/imputation_rulegraph.svg)

## Tools

The rules drive standard population-genetics tooling:

| Step | Tool |
|------|------|
| VCF manipulation, reference fixing | `bcftools` (incl. `+fixref`, `+liftover`) |
| QC, ancestry, format conversion | `plink` / `plink2`, `vcftools` |
| Liftover hg19 → GRCh38 | `CrossMap` |
| Strand/allele harmonisation | `GenotypeHarmonizer` |
| Phasing | `Eagle` |
| Imputation | `Minimac4` |

> **Tools are public and version-pinned.** The former monolithic image
> (`SNP_imputation_1000g_hg38.sif`) has been retired in favour of per-rule conda
> environments (`workflow/envs/`), each pinning a public, versioned tool; a few
> rules use pinned public biocontainers (e.g. `bcftools:1.21`). The two tools with
> no bioconda package (GenotypeHarmonizer, Minimac4 1.0.2) are fetched from public
> GitHub releases via `dvc import-url`. See issue #23.

## Reference data and provenance

All reference data is **public** and tracked via `dvc import-url` from public
sources (see the `.dvc` files under `resources/`), so it can be fetched without
access to any private registry:

- `resources/eQTLGenImpRef.tar.gz` — sceQTL-Gen / Powell Lab bundle (Minimac4
  imputation reference, 1000 Genomes 30x-GRCh38 phasing/panel, GRCh38 QC FASTA)
- `resources/1000G.tar.gz` — 1000 Genomes phase-3 plink (ancestry QC)
- `resources/bed/hg38exonsUCSC.bed` — hg38 exon BED
- `resources/liftover/GRCh37_to_GRCh38.chain.gz` — Ensembl GRCh37→GRCh38 liftover chain
- `resources/tools/` — GenotypeHarmonizer and Minimac4 (tools with no bioconda package)
- `resources/genomes/chr_map/`, `resources/genomes/GRCh38_chr.fai` — small, git-tracked

Fetch, then extract into the paths the config/rules expect:

```bash
# If you have access to the project's DVC remote:
dvc pull
# Or, without remote access, re-download straight from the public source URLs:
dvc update resources/eQTLGenImpRef.tar.gz.dvc resources/1000G.tar.gz.dvc \
           resources/bed/hg38exonsUCSC.bed.dvc \
           resources/liftover/GRCh37_to_GRCh38.chain.gz.dvc \
           resources/tools/GenotypeHarmonizer-1.4.23-dist.tar.gz.dvc \
           resources/tools/minimac4-1.0.2-Linux.sh.dvc

./prep.sh    # extract the archives into place (~160 GB); tarballs can then be deleted
```

This workflow was **adapted from the Powell Lab / sceQTL-Gen consortium
imputation pipeline**
([powellgenomicslab/SNP_imputation_1000g_hg38](https://github.com/powellgenomicslab/SNP_imputation_1000g_hg38)),
itself developed for the sceQTL-Gen WG1 pipeline. If you use this workflow,
please cite:

- van der Wijst *et al.* (2020), *eLife* — the sceQTL-Gen / single-cell eQTL
  reference approach.
- Neavin *et al.* (2024), *Genome Biology* — **Demuxafy**, the downstream
  demultiplexing framework these SNP profiles feed.
- The reference panels and tools above per their own citation requirements
  (1000 Genomes Project; Minimac4; Eagle; plink; bcftools; GenotypeHarmonizer;
  CrossMap).

See `CITATION.cff` for how to cite this workflow itself.

## Running

Each stage has a driver script that runs from the top of a super-project that
mounts this repo as a module (see
[brca_mega_atlas-data](https://github.com/swarbricklab/brca_mega_atlas-data)):

```bash
./modules/imputation/run_qc.sh           # QC stage
./modules/imputation/run_imputation.sh   # imputation stage
```

On NCI the driver uses the shared PBS profile carried as a nested submodule
(`profiles/global` → `swarbricklab/snakemake_config`). Check the repo out
recursively:

```bash
git submodule update --init --recursive
```

## Status

This repository is being prepared for public release alongside the manuscript.
See `docs/PUBLICATION_READINESS.md` for the roadmap and the open issues it maps
to.
