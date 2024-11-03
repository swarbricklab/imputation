import pysam
from pathlib import Path
import sys

# Redirect print statements to the Snakemake log file
sys.stdout = open(snakemake.log[0], 'w')

input_vcf = Path(snakemake.input['vcf'])
output_vcf = Path(snakemake.output['vcf'])

print('Opening the input VCF')
with pysam.VariantFile(input_vcf, 'r') as vcf_in:
    print('Creating a new VCF file with the same header')
    with pysam.VariantFile(output_vcf, 'w', header=vcf_in.header) as vcf_out:
        print('Modifying contig names and processing records')
        for rec in vcf_in.fetch():
            if rec.contig == '24':
                print('Renaming chromosome 24 to X')
                rec.contig = 'X'
            elif rec.contig.isdigit():
                print('Renaming chromosomes to add chr prefix if needed (e.g., 1 to chr1)')
                rec.contig = f'chr{rec.contig}'
            elif rec.contig == 'MT':
                print('Renaming chromosome MT to chrM')
                rec.contig = 'chrM'
            print('Writing modified record')
            vcf_out.write(rec)

sys.stdout.close()
