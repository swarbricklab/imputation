#!/bin/bash
# Touch empty placeholder inputs so `snakemake -n` can build both DAGs (QC +
# imputation) in CI without any real data. Mirrors config/test_{qc,imputation}.yaml.
set -eu
mk() { mkdir -p "$(dirname "$1")"; : > "$1"; }

# --- QC leaf inputs ---
mk .tests/integration/data/test.hg19.vcf.gz
mk .tests/integration/data/test.psam
mk .tests/integration/data/id_map.tsv
mk resources/1000g/all_phase3_filtered.pgen
mk resources/1000g/all_phase3_filtered.pvar
mk resources/1000g/all_phase3_filtered.psam

# --- imputation leaf inputs ---
mk .tests/integration/data/test.hg38.vcf.gz
# post_qc.psam is read at PARSE time for the Provided_Ancestry column -> real content
mkdir -p .tests/integration/results/post_qc
printf '#FID\tIID\tPAT\tMAT\tSEX\tProvided_Ancestry\n0\tS1\t0\t0\t1\tEUR\n0\tS2\t0\t0\t2\tAFR\n' \
    > .tests/integration/results/post_qc/post_qc.psam
mk .tests/integration/results/post_qc/post_qc.pgen
mk .tests/integration/results/post_qc/post_qc.pvar
mk resources/genomes/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted.vcf.gz
mk resources/genomes/hg38/ref_panel_QC/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi
mk resources/genomes/hg38/ref_genome_QC/Homo_sapiens.GRCh38.dna.primary_assembly.fa
mk resources/genomes/GRCh38_chr.fai
mk resources/genomes/hg38/phasing/genetic_map/genetic_map_hg38_withX.txt.gz
mk resources/genomes/chr_map/int2chr.txt
mk resources/genomes/chr_map/chr2int.txt
mk resources/bed/hg38exonsUCSC.bed
mk resources/liftover/GRCh37_to_GRCh38.chain.gz
for c in $(seq 1 22); do
  mk "resources/genomes/hg38/phasing/phasing_reference/chr${c}.bcf"
  mk "resources/reference/hg38/imputation/chr${c}.m3vcf.gz"
done
echo "placeholders touched"
