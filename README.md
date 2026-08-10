# Imputation

This workflow starts with a VCF file containing measured genotypes anmd imputes missing variants based on data from the 1000 Genomes Project.
Strictly speaking, this repo contains not one but two workflows: 1) a quality control workflow that filters SNPs and checks sex and ancestry annotations against the provided genotype information, and 2) an imputation workflow that imputes missing variants using data from the 1000 Genomes Project as a reference.
Imputed variants are further filtered to exclude indels and SNPs outside exonic regions, as well as variants with low minor allele frequencies (MAF) or poor correlation scores.

These workflows are intermediate steps in SNP demultiplexing with demuxafy, which uses the SNP profiles in the final VCF file to assign cells to donors in multiplexed pool of cells that have been captured together with 10X Chromium.

# Overview

The `plink_qc` and `imputation` workflows in this repo are the second and third stages in the `demuxafy` SNP demux process, as shown below.

TODO: insert DAG


Rule graphs for each stage are shown below. 
See the individual rule definitions to understand the function of each rule.

## Plink QC

This stage starts with the following outputs from the [genotyping](url):
- VCF: genotype calls for all samples, based on raw SNP microarray data
- psam: a sample description file, ready for processing by `plink`
- id map: `plink` is fussy about sample ids, so this file maps the original sample ids to the adjusted ids in the `.psam` file

The output of this stage is an updated set of plink files:
- pgen: genotype information in plink format
- pvar: ?
- psam: an updated version of the sample description file, in which sex annoations have been double checked and ancestries have been assigned based on the 1000 Genomes reference
- unrelated_sample_ids.txt: a maximal unrelated subset of the post-QC samples, selected by KING kinship (`--king-cutoff`)
- derived_sample_metadata.csv: per-sample sex, inferred ancestry, kinship group and genotyping call rate

### Related samples

Samples that fail the sex or ancestry checks are not removed; their annotations are corrected instead (sex to the genotype-inferred value, ancestry to the PCA assignment). Related samples are also **not** removed: the final VCF exists for SNP demultiplexing, so every array a donor has must stay available for selection downstream.
Instead, the KING keep-list (`unrelated_sample_ids.txt`) defines the unrelated subset on which the imputation stage computes its frequency-based variant statistics (HWE, MAF), so duplicate samples cannot distort them; the resulting variant filters are then applied to all samples (see rule `filter_preimpute_vcf`). Per-variant missingness is deliberately evaluated over all samples, since every array is an independent observation of probe quality.
In the metadata, samples of the same donor (kinship at or above `king_duplicate_cutoff`, i.e. repeat arrays or identical twins) share a `kinship_group`, and `genotype_call_rate` gives a basis for choosing among them.
The default cutoff (0.3) aims to accomodate a loss of heterozygosity depressesing the kinship estimate between a donor's tumour and blood arrays, while staying well above the first-degree relative expectation of 0.25.

The rule graph is as follows:

![QC Rulegraph](docs/qc_rulegraph.svg)

## Imputation

![Imputation Rulegraph](docs/imputation_rulegraph.svg)
