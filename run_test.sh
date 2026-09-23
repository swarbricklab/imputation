
#! /bin/bash

set -e
eval  "$(conda shell.bash hook)"
conda activate snakemake_7.32.4

module=imputation

if [[ "$(hostname)" == *"nci"* ]]; then
    echo "Running on NCI"
    global_profile="--profile profiles/global/nci"
    workflow_profile="--workflow-profile profiles/workflow"
    module load singularity
    mkdir -p logs/joblogs
else
    echo "WARNING: Unknown host: $(hostname)"
    echo "WARNING: No known global profile for this host"
    echo "WARNING: Running without global profile"
    echo "See https://github.com/swarbricklab/snakemake_config/blob/main/README.md"
    global_profile=""
    workflow_profile="--workflow-profile profiles/workflow"
fi

# QC
snakemake $global_profile $workflow_profile \
    --snakefile workflow/Snakefile_qc \
    --configfile config/test_qc.yaml \
    $@

# Imputation
snakemake $global_profile $workflow_profile \
    --snakefile workflow/Snakefile_imputation \
    --configfile config/test_imputation.yaml \
    $@
