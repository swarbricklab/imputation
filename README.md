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

The rule graph is as follows:

TODO: insert rule graph 


