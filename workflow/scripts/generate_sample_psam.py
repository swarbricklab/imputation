import pysam
from pathlib import Path
import pandas as pd
import sys

# Redirect print statements to the Snakemake log file
sys.stdout = open(snakemake.log[0], 'w')

# Get input/output paths
input_vcf = Path(snakemake.input['vcf'])
psam_output = Path(snakemake.output['psam'])

print('Extracting sample IDs (IIDs) from the VCF header')
with pysam.VariantFile(input_vcf, 'r') as vcf_in:
    sample_names = list(vcf_in.header.samples)

# Create the metadata.psam file with default values for FID, PAT, MAT, SEX, and Provided_Ancestry
print('Creating the metadata.psam DataFrame with default values')
psam_df = pd.DataFrame({
    '#FID': [0] * len(sample_names),
    'IID': sample_names,
    'PAT': [0] * len(sample_names),
    'MAT': [0] * len(sample_names),
    'SEX': [0] * len(sample_names),
    'Provided_Ancestry': ['NONE'] * len(sample_names)
})
print('Writing metadata.psam')
psam_df.to_csv(psam_output, sep='\t', index=False)

sys.stdout.close()