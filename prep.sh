#! /bin/bash

# Prepare / check the imputation reference data and the two non-bioconda tools.
#
# The reference archives and tools are tracked publicly via `dvc import-url` (see
# the README and resources/**/*.dvc). This script unpacks them into the paths the
# workflow config / rules expect. Run it, after `dvc pull`, from the top of the
# repo (standalone) or of the super-project that mounts imputation as a module --
# the same place you run run_qc.sh / run_imputation.sh.
#
#   resources/eQTLGenImpRef.tar.gz  (sceQTL-Gen / Powell Lab bundle, ~36 GB) ->
#       resources/genomes/hg38/{phasing,ref_genome_QC,ref_panel_QC}
#       resources/reference/hg38/imputation           (Minimac4 reference)
#   resources/1000G.tar.gz          (1000 Genomes phase-3 plink, ancestry QC) ->
#       resources/1000g/all_phase3_filtered.{pgen,pvar,psam}
#   resources/tools/GenotypeHarmonizer-1.4.23-dist.tar.gz ->
#       resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar (run on openjdk)
#   resources/tools/minimac4-1.0.2-Linux.sh -> resources/tools/minimac4  (static binary;
#       bioconda only ships Minimac4 4.x, which needs .msav not our 1.x .m3vcf panel)
#
# Usage:
#   dvc pull                                          # fetch resources/**/*.dvc
#   ./prep.sh                                         # unpack (default: --what all)
#   ./prep.sh --what check --configfile config/test_qc.yaml   # read-only preflight
#
# Options:
#   --what all|check   all: unpack the archives (default). check: a read-only
#                      preflight that verifies the references, tools and the
#                      config's inputs are present, and exits non-zero if not, so
#                      it can gate a run. Fetches and unpacks nothing.
#   --configfile PATH  config to read the deps/refs paths from (required for check).
#   --force            re-unpack even if the outputs already look complete.
#
# Small references are git-tracked and need no preparation:
#   resources/genomes/chr_map/{int2chr,chr2int}.txt   (chromosome-name maps)
#   resources/genomes/GRCh38_chr.fai                  (chr-named GRCh38 .fai)

set -euo pipefail

what="all"
configfile=""
force="false"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --what)       what="${2:-}"; shift 2 ;;
        --configfile) configfile="${2:-}"; shift 2 ;;
        --force)      force="true"; shift ;;
        -h|--help)    sed -n '3,45p' "$0" | cut -c 3-; exit 0 ;;
        *)            echo "ERROR: unknown argument: $1" >&2; exit 1 ;;
    esac
done
case "$what" in all|check) ;; *) echo "ERROR: --what must be all or check (got: $what)" >&2; exit 1 ;; esac

die() { echo "ERROR: $*" >&2; exit 1; }

# Make sure python can read YAML (only needed for the check); fall back to the
# snakemake env, which has PyYAML.
ensure_pyyaml() {
    python3 -c "import yaml" 2>/dev/null && return 0
    eval "$(conda shell.bash hook)"
    conda activate snakemake_7.32.4
}


# Read-only preflight: confirm the references, tools and the config's inputs are
# present before a run, without fetching or unpacking anything. Reuses the same
# config so it checks exactly what the workflow will look for, and exits non-zero
# if anything is missing so a super-project can gate a run on it.
preflight() {
    [[ -n "$configfile" ]] || die "give --configfile (see --help)."
    [[ -f "$configfile" ]] || die "config file not found: $configfile"
    ensure_pyyaml

    local fails=0
    ok()  { echo "  ok    $1"; }
    bad() { echo "  FAIL  $1"; fails=$((fails + 1)); }

    echo "Preflight for $configfile (read-only; fetches and unpacks nothing)"
    echo

    # 1. conda/mamba: the rules run in per-rule conda envs.
    echo "Conda (rules run in workflow/envs/ conda environments):"
    if command -v mamba > /dev/null; then ok "mamba available"
    elif command -v conda > /dev/null; then ok "conda available"
    else bad "neither mamba nor conda found -- needed to build the rule envs"; fi
    echo

    # 2. Inputs and references: every string path under deps: and refs:.
    echo "Inputs and references (deps + refs in $configfile):"
    local paths
    paths="$(python3 - "$configfile" <<'PY'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1])) or {}
def leaves(x):
    if isinstance(x, dict):
        for v in x.values(): yield from leaves(v)
    elif isinstance(x, list):
        for v in x: yield from leaves(v)
    elif isinstance(x, str):
        yield x
for sec in ("deps", "refs"):
    for p in leaves(d.get(sec) or {}):
        print(p)
PY
)"
    local p
    while IFS= read -r p; do
        [[ -n "$p" ]] || continue
        if [[ -e "$p" ]]; then ok "$p"
        else bad "missing: $p -- 'dvc pull' + './prep.sh' (unpack), or check the path"; fi
    done <<< "$paths"
    echo

    # 3. Non-bioconda tools (only relevant to the imputation stage, which sets refs.impute).
    if python3 - "$configfile" <<'PY'
import sys, yaml
d = yaml.safe_load(open(sys.argv[1])) or {}
sys.exit(0 if (d.get("refs") or {}).get("impute") else 1)
PY
    then
        echo "Tools (GenotypeHarmonizer, Minimac4 -- unpacked by prep.sh):"
        local t
        for t in resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar resources/tools/minimac4; do
            [[ -s "$t" ]] && ok "$t" || bad "missing: $t -- run './prep.sh' to unpack"
        done
        echo
    fi

    if [[ "$fails" -eq 0 ]]; then
        echo "Preflight OK: all inputs, references and tools are in place."
    else
        die "preflight found $fails problem(s) above; resolve them before running."
    fi
}


if [[ "$what" == "check" ]]; then
    preflight
    exit 0
fi

# ============================ unpack (--what all) ============================

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
    bash "$mm" --skip-license --prefix="$inst" > /dev/null
    bin="$(find "$inst" -type f -name minimac4 | head -1)"
    [[ -n "$bin" ]] || { echo "ERROR: minimac4 binary not found after installing $mm" >&2; exit 1; }
    cp "$bin" resources/tools/minimac4
    chmod +x resources/tools/minimac4
    rm -rf "$inst"; trap - EXIT
    echo "  -> resources/tools/minimac4"
fi

echo "Done."
