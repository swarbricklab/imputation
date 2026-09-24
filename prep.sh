#! /bin/bash

# Prepare the imputation reference data from its public source.
#
# The bulk of the reference set -- the Minimac4 imputation reference, the 1000
# Genomes 30x-GRCh38 phasing/panel data, and the GRCh38 QC FASTA -- is the
# sceQTL-Gen / Powell Lab bundle, tracked publicly via `dvc import-url` as
# resources/eQTLGenImpRef.tar.gz (see resources/eQTLGenImpRef.tar.gz.dvc and the
# README). This script extracts that bundle into the paths the workflow config
# expects:
#
#   resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}
#   resources/reference/hg38/imputation
#
# Run it once, after `dvc pull resources/eQTLGenImpRef.tar.gz.dvc` (or the
# `dvc import-url`) has fetched the ~36 GB bundle. It expands to ~120 GB.
#
# Usage:
#   ./prep.sh [--force]
#
# Not handled here yet -- still separate references (see the README and issue
# #22): resources/1000g (1000G phase-3 plink for ancestry QC), resources/bed
# (hg38 exon BED), resources/genomes/chr_map.

set -euo pipefail

force="false"
[[ "${1:-}" == "--force" ]] && force="true"

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$here"

bundle="resources/eQTLGenImpRef.tar.gz"
if [[ ! -s "$bundle" ]]; then
    echo "ERROR: $bundle is not present." >&2
    echo "Fetch it first, e.g.:  dvc pull resources/eQTLGenImpRef.tar.gz.dvc" >&2
    exit 1
fi

# Already extracted? (unless --force)
if [[ "$force" != "true" \
      && -d resources/reference/hg38/imputation \
      && -d resources/genomes/hg38/ref_panel_QC \
      && -d resources/genomes/hg38/phasing \
      && -d resources/genomes/hg38/ref_genome_QC ]]; then
    echo "Imputation reference already extracted. Nothing to do (use --force)."
    exit 0
fi

echo "Extracting $bundle (expands to ~120 GB) ..."
tmp="$(mktemp -d resources/.eqtlgen.XXXXXX)"
trap 'rm -rf "$tmp"' EXIT
tar xzf "$bundle" -C "$tmp"

[[ -d "$tmp/hg38" ]] \
    || { echo "ERROR: unexpected bundle layout (expected hg38/ at the root)." >&2; exit 1; }

mkdir -p resources/genomes/hg38 resources/reference/hg38

# phasing / ref_genome_QC / ref_panel_QC -> resources/genomes/hg38/
for d in phasing ref_genome_QC ref_panel_QC; do
    [[ -e "$tmp/hg38/$d" ]] || { echo "ERROR: bundle missing hg38/$d" >&2; exit 1; }
    rm -rf "resources/genomes/hg38/$d"
    mv "$tmp/hg38/$d" "resources/genomes/hg38/$d"
done

# imputation (Minimac4 reference) -> resources/reference/hg38/imputation
[[ -e "$tmp/hg38/imputation" ]] || { echo "ERROR: bundle missing hg38/imputation" >&2; exit 1; }
rm -rf resources/reference/hg38/imputation
mv "$tmp/hg38/imputation" resources/reference/hg38/imputation

echo "Done. Imputation reference extracted:"
echo "  resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}"
echo "  resources/reference/hg38/imputation"
