
#!/usr/bin/env python
shell.executable('bash')



# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = imputation/"subset_ancestry/{ancestry}_subset.pgen",
        psam = imputation/"subset_ancestry/{ancestry}_subset.psam",
        pvar = imputation/"subset_ancestry/{ancestry}_subset.pvar"
    output:
        bed = imputation/"crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = imputation/"crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = imputation/"crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = imputation/"crossmapped/{ancestry}_crossmap_input.bed",
        outbed = imputation/"crossmapped/{ancestry}_crossmap_output.bed",
        excluded_ids = imputation/"crossmapped/{ancestry}_excluded_ids.txt"
    params:
        in_plink = imputation/"subset_ancestry/{ancestry}_subset",
        out = imputation/"crossmapped/{ancestry}_crossmapped_plink",
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}
        plink2 --pfile {params.in_plink} --exclude {output.excluded_ids} --make-bed --output-chr MT --out {params.out}
        awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        pgen = imputation/"crossmapped/{ancestry}_crossmapped_plink.bed",
        psam = imputation/"crossmapped/{ancestry}_crossmapped_plink.bim",
        pvar = imputation/"crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    params:
        infile = imputation/"crossmapped/{ancestry}_crossmapped_plink",
        out = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --bfile {params.infile} --make-bed --max-alleles 2 --output-chr MT --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi"
    output:
        bed = imputation/"harmonize_hg38/{ancestry}.bed",
        bim = imputation/"harmonize_hg38/{ancestry}.bim",
        fam = imputation/"harmonize_hg38/{ancestry}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
    params:
        infile = imputation/"crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = imputation/"harmonize_hg38/{ancestry}",
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        java -Xmx{resources.java_mem}g -jar {params.jar}\
            --input {params.infile}\
            --inputType PLINK_BED\
            --ref {input.vcf}\
            --refType VCF\
            --update-id\
            --output {params.out}
        """


rule plink_to_vcf:
    input:
        bed = imputation/"harmonize_hg38/{ancestry}.bed",
        bim = imputation/"harmonize_hg38/{ancestry}.bim",
        fam = imputation/"harmonize_hg38/{ancestry}.fam"
    output:
        data_vcf_gz = imputation/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = imputation/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    params:
        infile = imputation/"harmonize_hg38/{ancestry}",
        out = imputation/"harmonize_hg38/{ancestry}_harmonised_hg38"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --bfile {params.infile} --recode vcf id-paste=iid --chr 1-22 --out {params.out}

        bgzip {params.out}.vcf
        bcftools index {output.data_vcf_gz}
        """


rule vcf_fixref_hg38:
    input:
        fasta = config['ref']['dir'] + config['ref']['fasta'],
        vcf = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi",
        data_vcf = imputation/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf = imputation/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = imputation/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} | \
        bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}

        #Index
        bcftools index {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = imputation/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf = imputation/"filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = imputation/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = imputation/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    params:
        maf = lambda wildcards: float(maf_df["MAF"][maf_df.Ancestry == wildcards.ancestry].values),
        missing = config["imputation"]["snp_missing_pct"],
        hwe = config["imputation"]["snp_hwe"]
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        #Add tags
        export BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins
        bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} |\
        bcftools filter -e 'REF="N" | REF="I" | REF="D"' |\
        bcftools filter -e "ALT='.'" |\
        bcftools norm -d all |\
        bcftools norm -m+any |\
        bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        #Index the output file
        bcftools index {output.filtered_vcf}
        """

rule het:
    input:
        vcf = imputation/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf = temp(imputation/"het/{ancestry}_filtered_temp.vcf"),
        inds = imputation/"het/{ancestry}_het_failed.inds",
        het = imputation/"het/{ancestry}_het.het",
        passed = imputation/"het/{ancestry}_het_passed.inds",
        passed_list = imputation/"het/{ancestry}_het_passed.txt"
    params:
        het_base = imputation/"het/{ancestry}_het",
        script = "/opt/SNP_imputation_1000g_hg38/Imputation/scripts/filter_het.R",
        hwe = imputation/"hwe/{ancestry}_hwe",
        out = imputation/"het/{ancestry}_het"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        gunzip -c {input.vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}
        vcftools --vcf {output.tmp_vcf} --het --out {params.het_base}
        Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list = imputation/"het/{ancestry}_het_passed.txt",
        vcf = imputation/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz",
        index = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"]
    params:
        hwe = imputation/"hwe/{ancestry}_hwe",
        out = imputation/"het_filter/{ancestry}_het_filter"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}

        #Index the output file
        bcftools index {output.vcf}
        """


rule calculate_missingness:
    input:
        filtered_vcf = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(imputation/"filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss = imputation/"filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals = imputation/"genotype_donor_annotation/{ancestry}_individuals.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"]
    params:
        out = imputation/"filter_preimpute_vcf/{ancestry}_genotypes"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        gunzip -c {input.filtered_vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}

        vcftools --gzvcf {output.tmp_vcf} --missing-indv --out {params.out}

        bcftools query -l {input.filtered_vcf} >> {output.individuals}
        """


rule split_by_chr:
    input:
        filtered_vcf = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = imputation/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = imputation/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        index = imputation/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """


rule eagle_prephasing:
    input:
        vcf = imputation/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = Path(config['ref']['dir'])/config['ref']['genetic_map'],
        phasing_file = Path(config['ref']['dir'])/config['ref']['phasing_dir']/"chr{chr}.bcf"
    output:
        vcf = imputation/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    params:
        out = imputation/"eagle_prephasing/{ancestry}_chr{chr}_phased"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix={params.out} \
            --numThreads={threads}
        """


rule minimac_imputation:
    input:
        vcf = imputation/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file = config['ref']['dir'] + config['ref']['impute_dir'] + "/chr{chr}.m3vcf.gz"
    output:
        imputation/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    params:
        out = imputation/"minimac_imputed/{ancestry}_chr{chr}",
        minimac4 = "/opt/bin/minimac4",
        chunk_length = config["imputation"]["chunk_length"]
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        {params.minimac4} --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix {params.out} \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads} \
            --ChunkLengthMb {params.chunk_length}
        """


rule combine_vcfs_ancestry:
    input:
        vcfs = lambda wildcards: expand(imputation/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", chr = chromosomes, ancestry = ancestry_subsets)
    output:
        combined = imputation/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        ind = imputation/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    params:
        files_begin = imputation/"minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        bcftools concat -Oz {params.files_begin} > {output.combined}
        bcftools index {output.combined}
        """


rule combine_vcfs_all:
    input:
        vcfs = lambda wildcards: expand(imputation/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets)
    output:
        combined = imputation/"vcf_all_merged/imputed_hg38.vcf.gz",
        ind = imputation/"vcf_all_merged/imputed_hg38.vcf.gz.csi"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        if [[ $(ls -l {input.vcfs} | wc -l) > 1 ]]
        then
            bcftools merge -Oz {input.vcfs} > {output.combined}
        else
            cp {input.vcfs} {output.combined}
        fi
        bcftools index {output.combined}
        """
