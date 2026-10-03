
#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

module=imputation

if [[ "$(hostname)" == *"nci"* ]]; then
    echo "Running on NCI"
    global_profile="--profile modules/$module/profiles/global/nci"
    workflow_profile="--workflow-profile modules/$module/profiles/workflow "
    module load singularity
    mkdir -p logs/joblogs
else
    echo "WARNING: Unknown host: $(hostname)"
    echo "WARNING: No known global profile for this host"
    echo "WARNING: Running without global profile"
    echo "See https://github.com/swarbricklab/snakemake_config/blob/main/README.md"
    global_profile=""
fi

# Preflight: verify references, tools and inputs are present before submitting
# jobs. Set SKIP_PREFLIGHT=1 to bypass.
if [[ -z "${SKIP_PREFLIGHT:-}" ]]; then
    modules/$module/prep.sh --what check --configfile config/$module/config_qc.yaml
else
    echo "SKIP_PREFLIGHT set -- skipping the dependency preflight."
fi

snakemake $global_profile $workflow_profile \
    --snakefile modules/$module/workflow/Snakefile_qc \
    --configfile config/$module/config_qc.yaml \
    $@
