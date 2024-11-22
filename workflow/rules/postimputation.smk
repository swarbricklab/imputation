# =============

# ============= 

# ==================================
# vcf_filter
# ----------------------------------
# Preparing the VCF file for downstream SNP demux tools
# ------------
# TODO: 
#   - more detailed description of what this rule does
#   - break-down into more discrete steps (the shell commands do a lot of different filtering)
# ==================================
rule vcf_filter:
    input:
        vcf=imputation/"vcf_all_merged/imputed_hg38.vcf.gz"
    output: 
        # TODO: should the output directory be "postimputation" rather than imputation?
        vcf_MAF_R2_filtered=imputation/"demuxafy_input/Merged_Imputed_R2_0.3_MAF0.05.vcf",
        vcf_merged_recode=imputation/"demuxafy_input/Merged_Imputed_R2_0.3_MAF0.05_Exon.vcf.recode.vcf",
        vcf_merged_recode_X=imputation/"demuxafy_input/Merged_Imputed_R2_0.3_MAF0.05_Exon_chrX.vcf",
        vcf_merged_recode_X_24=imputation/"demuxafy_input/Merged_Imputed_R2_0.3_MAF0.05_Exon_chrX_24.vcf",
        vcf_reheader=imputation/"demuxafy_input/reheader_Merged_Imputed_R2_0.3_MAF0.05_Exon.vcf"
    params:
        out_dir=imputation/"demuxafy_input/",
        BED=config["exons_bed"],
        FAI=config["genome_index"]
    log:
        logs/"imputation/demuxafy_input/vcf_filter.log"
    container:
        config["containers"]["SNP_imputation"] 
    shell:
        """
        # filter MAF and R2 values
        cmd="bcftools filter --include 'MAF>=0.05 & R2>=0.3' -O v --output {output.vcf_MAF_R2_filtered} {input.vcf}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        # filter exon 
        # NOTE/TODO: the --out filename assignment is a bit of a hack due to the --recode param adding a .recode.vcf file suffix
        cmd="vcftools \
        --vcf {output.vcf_MAF_R2_filtered} \
        --max-alleles 2 \
        --remove-indels \
        --bed {params.BED} \
        --recode \
        --recode-INFO-all \
        --out {params.out_dir}Merged_Imputed_R2_0.3_MAF0.05_Exon.vcf" 
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        cmd="sed -i '/##contig=<ID=/d' {output.vcf_merged_recode}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        cmd="sed 's/^24/X/g' {output.vcf_merged_recode}  > {output.vcf_merged_recode_X}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        cmd="grep -v 2147483648 {output.vcf_merged_recode_X} > {output.vcf_merged_recode_X_24}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        cmd="sed -i '/UNKNOWNPOSITION/d' {output.vcf_merged_recode_X_24}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        # #### TO match the BAM and VCF (they should be in same order)
        # TODO: Is this a correct description? Where is the bam file here?
        cmd="bcftools reheader -f {params.FAI} {output.vcf_merged_recode_X_24} > {output.vcf_reheader}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        echo "vcf_filter COMPLETE" >> {log}
        """

# ==================================
# vcf_filter2
# ----------------------------------
# Second step to preparing the VCF file for downstream SNP demux tools
#  - Sorts the imputed known-genotypes VCF file
#  - Outputs both an int and chr encoded version of the VCF file
# ==================================
rule vcf_filter2:
    input:
        vcf=imputation/"demuxafy_input/reheader_Merged_Imputed_R2_0.3_MAF0.05_Exon.vcf"
    output: 
        vcf_int_encoded=imputation/"demuxafy_input/reheader_Merged_Imputed_R2_0.3_MAF0.05_Exon_sorted.int_encoded.vcf",
        vcf_chr_encoded=imputation/"demuxafy_input/reheader_Merged_Imputed_R2_0.3_MAF0.05_Exon_sorted.chr_encoded.vcf"
    log:
        logs/"imputation/demuxafy_input/vcf_filter2.log"
    container:
        config["containers"]["SNP_imputation"] 
    shell:
        """
        cmd="java -jar /opt/picard/build/libs/picard.jar SortVcf \
        I={input.vcf} \
        O={output.vcf_int_encoded}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1

        cmd="sed -e '/^#/! s/^\([0-9XY]\)/chr\\1/' \
            -e '/^##contig=<ID=/ s/\(##contig=<ID=\)\(.*\)/\\1chr\\2/' \
            {output.vcf_int_encoded} \
            > {output.vcf_chr_encoded}"
        echo $cmd >> {log}
        eval $cmd >> {log} 2>&1
        """       

# =============
# ready to run demuxafy tools
# ============= 