#! /bin/bash

# Prepare the imputation reference data and the two non-bioconda tools from
# their public sources.
#
# Everything below is tracked publicly via `dvc import-url` (see the README and
# resources/**/*.dvc). This script unpacks the archives into the paths the
# workflow config / rules expect. Run it once, after `dvc pull`.
#
#   resources/eQTLGenImpRef.tar.gz  (sceQTL-Gen / Powell Lab bundle, ~36 GB) ->
#       resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}
#       resources/reference/hg38/imputation           (Minimac4 reference)
#
#   resources/1000G.tar.gz          (1000 Genomes phase-3 plink, ancestry QC) ->
#       resources/1000g/all_phase3_filtered.{pgen,pvar,psam}
#
# Tools not on bioconda (every other tool comes from a conda env, workflow/envs/):
#   resources/tools/GenotypeHarmonizer-1.4.23-dist.tar.gz ->
#       resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar  (run on openjdk)
#   resources/tools/minimac4-1.0.2-Linux.sh -> resources/tools/minimac4   (static binary;
#       bioconda only ships Minimac4 4.x, which needs .msav not our 1.x .m3vcf panel)
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

# All paths below are relative to the current directory, so run this from the top
# of the repo (standalone) or of the super-project that mounts it as a module --
# the same place you run run_qc.sh / run_imputation.sh.

# --- eQTLGen imputation reference bundle -----------------------------------
bundle="resources/eQTLGenImpRef.tar.gz"
if [[ "$force" != "true" \
      && -d resources/reference/hg38/imputation \
      && -d resources/genomes/hg38/ref_panel_QC \
      && -d resources/genomes/hg38/phasing \
      && -d resources/genomes/hg38/ref_genome_QC ]]; then
    echo "eQTLGen reference already extracted."
else
    # Only require the tarball when extraction is actually needed, so the large
    # tarballs can be deleted once the data is extracted.
    if [[ ! -s "$bundle" ]]; then
        echo "ERROR: $bundle is not present. Fetch it first: dvc pull $bundle.dvc (or dvc update $bundle.dvc)" >&2
        exit 1
    fi
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
if [[ "$force" != "true" \
      && -s resources/1000g/all_phase3_filtered.pgen \
      && -s resources/1000g/all_phase3_filtered.pvar \
      && -s resources/1000g/all_phase3_filtered.psam ]]; then
    echo "1000G reference already extracted."
else
    if [[ ! -s "$kg" ]]; then
        echo "ERROR: $kg is not present. Fetch it first: dvc pull $kg.dvc (or dvc update $kg.dvc)" >&2
        exit 1
    fi
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

# --- Tools that are not on bioconda (provisioned from public releases) ------
# GenotypeHarmonizer (Java jar) and Minimac4 1.0.2 (static binary) are tracked
# via dvc import-url (resources/tools/*.dvc) and unpacked here. Every other tool
# comes from a conda env (workflow/envs/); GenotypeHarmonizer runs on the
# openjdk env and Minimac4 is a static binary.
gh="resources/tools/GenotypeHarmonizer-1.4.23-dist.tar.gz"
if [[ ! -s "$gh" ]]; then
    echo "ERROR: $gh is not present. Fetch it first: dvc pull $gh.dvc" >&2
    exit 1
fi
if [[ "$force" == "true" || ! -s resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar ]]; then
    echo "Unpacking GenotypeHarmonizer ..."
    tar xzf "$gh" -C resources/tools
    [[ -s resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar ]] \
        || { echo "ERROR: GenotypeHarmonizer.jar not found after unpacking $gh" >&2; exit 1; }
    echo "  -> resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
fi

mm="resources/tools/minimac4-1.0.2-Linux.sh"
if [[ ! -s "$mm" ]]; then
    echo "ERROR: $mm is not present. Fetch it first: dvc pull $mm.dvc" >&2
    exit 1
fi
if [[ "$force" == "true" || ! -x resources/tools/minimac4 ]]; then
    echo "Installing Minimac4 1.0.2 ..."
    inst="$(mktemp -d resources/tools/.mm4.XXXXXX)"
    trap 'rm -rf "$inst"' EXIT
    # self-extracting CMake installer
    bash "$mm" --skip-license --prefix="$inst" > /dev/null
    bin="$(find "$inst" -type f -name minimac4 | head -1)"
    [[ -n "$bin" ]] || { echo "ERROR: minimac4 binary not found after installing $mm" >&2; exit 1; }
    cp "$bin" resources/tools/minimac4
    chmod +x resources/tools/minimac4
    rm -rf "$inst"; trap - EXIT
    echo "  -> resources/tools/minimac4"
fi

echo "Done."
