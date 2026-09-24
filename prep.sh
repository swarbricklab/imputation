#! /bin/bash

# Prepare the imputation reference data from its public sources.
#
# Two reference archives are tracked publicly via `dvc import-url` (see the
# README and resources/*.dvc). This script extracts them into the paths the
# workflow config expects. Run it once, after `dvc pull` has fetched the
# archives.
#
#   resources/eQTLGenImpRef.tar.gz  (sceQTL-Gen / Powell Lab bundle, ~36 GB) ->
#       resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}
#       resources/reference/hg38/imputation           (Minimac4 reference)
#
#   resources/1000G.tar.gz          (1000 Genomes phase-3 plink, ancestry QC) ->
#       resources/1000g/all_phase3_filtered.{pgen,pvar,psam}
#
# Usage:
#   dvc pull                 # fetch resources/*.tar.gz (and the exon BED)
#   ./prep.sh [--force]
#
# Small references are git-tracked and need no preparation:
#   resources/genomes/chr_map/{int2chr,chr2int}.txt   (chromosome-name maps)
#   resources/genomes/GRCh38_chr.fai                  (chr-named GRCh38 .fai,
#                                                      for bcftools reheader)
# The exon BED (resources/bed/hg38exonsUCSC.bed) is a single dvc import-url file
# and needs no extraction.

set -euo pipefail

force="false"
[[ "${1:-}" == "--force" ]] && force="true"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

# --- eQTLGen imputation reference bundle -----------------------------------
bundle="resources/eQTLGenImpRef.tar.gz"
if [[ ! -s "$bundle" ]]; then
    echo "ERROR: $bundle is not present. Fetch it first: dvc pull $bundle.dvc" >&2
    exit 1
fi
if [[ "$force" != "true" \
      && -d resources/reference/hg38/imputation \
      && -d resources/genomes/hg38/ref_panel_QC \
      && -d resources/genomes/hg38/phasing \
      && -d resources/genomes/hg38/ref_genome_QC ]]; then
    echo "eQTLGen reference already extracted."
else
    echo "Extracting $bundle (expands to ~120 GB) ..."
    tmp="$(mktemp -d resources/.eqtlgen.XXXXXX)"
    trap 'rm -rf "$tmp"' EXIT
    tar xzf "$bundle" -C "$tmp"
    [[ -d "$tmp/hg38" ]] || { echo "ERROR: unexpected bundle layout (no hg38/)." >&2; exit 1; }
    mkdir -p resources/genomes/hg38 resources/reference/hg38
    for d in phasing ref_genome_QC ref_panel_QC; do
        [[ -e "$tmp/hg38/$d" ]] || { echo "ERROR: bundle missing hg38/$d" >&2; exit 1; }
        rm -rf "resources/genomes/hg38/$d"; mv "$tmp/hg38/$d" "resources/genomes/hg38/$d"
    done
    [[ -e "$tmp/hg38/imputation" ]] || { echo "ERROR: bundle missing hg38/imputation" >&2; exit 1; }
    rm -rf resources/reference/hg38/imputation
    mv "$tmp/hg38/imputation" resources/reference/hg38/imputation
    rm -rf "$tmp"; trap - EXIT
    echo "  -> resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}"
    echo "  -> resources/reference/hg38/imputation"
fi

# --- 1000 Genomes phase-3 plink (ancestry QC) ------------------------------
kg="resources/1000G.tar.gz"
if [[ ! -s "$kg" ]]; then
    echo "ERROR: $kg is not present. Fetch it first: dvc pull $kg.dvc" >&2
    exit 1
fi
if [[ "$force" != "true" \
      && -s resources/1000g/all_phase3_filtered.pgen \
      && -s resources/1000g/all_phase3_filtered.pvar \
      && -s resources/1000g/all_phase3_filtered.psam ]]; then
    echo "1000G reference already extracted."
else
    echo "Extracting $kg ..."
    tmp="$(mktemp -d resources/.1000g.XXXXXX)"
    trap 'rm -rf "$tmp"' EXIT
    tar xzf "$kg" -C "$tmp"
    mkdir -p resources/1000g
    # Place the plink fileset wherever it sits inside the archive.
    moved=0
    for ext in pgen pvar psam; do
        f="$(find "$tmp" -type f -name "all_phase3_filtered.$ext" | head -1)"
        [[ -n "$f" ]] || { echo "ERROR: all_phase3_filtered.$ext not found in $kg" >&2; exit 1; }
        mv "$f" "resources/1000g/all_phase3_filtered.$ext"; moved=$((moved+1))
    done
    rm -rf "$tmp"; trap - EXIT
    echo "  -> resources/1000g/all_phase3_filtered.{pgen,pvar,psam} ($moved files)"
fi

echo "Done."
