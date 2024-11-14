rule subset_ancestry:
    input:
        psam = out_dir/"update_sex_ancestry/update_sex.psam"
    output:
        keep = out_dir/"subset_ancestry/{ancestry}_individuals.psam",
        pgen = out_dir/"subset_ancestry/{ancestry}_subset.pgen",
        psam = out_dir/"subset_ancestry/{ancestry}_subset.psam",
        pvar = out_dir/"subset_ancestry/{ancestry}_subset.pvar"
    params:
        infile = out_dir/"update_sex_ancestry/update_sex",
        out = out_dir/"subset_ancestry/{ancestry}_subset"
    container:
        config['deps']['container']
    shell:
        """
        grep {wildcards.ancestry} {input.psam} > {output.keep}
        plink2 --threads {threads} \
            --pfile {params.infile} \
            --keep {output.keep}  \
            --max-alleles 2 \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out {params.out}
        """

# Converts BIM to BED and converts the BED file via CrossMap.
# Finds excluded SNPs and removes them from the original plink file.
# Then replaces the BIM with CrossMap's output.
rule crossmap:
    input:
        pgen = out_dir/"subset_ancestry/{ancestry}_subset.pgen",
        psam = out_dir/"subset_ancestry/{ancestry}_subset.psam",
        pvar = out_dir/"subset_ancestry/{ancestry}_subset.pvar"
    output:
        bed = out_dir/"crossmapped/{ancestry}_crossmapped_plink.bed",
        bim = out_dir/"crossmapped/{ancestry}_crossmapped_plink.bim",
        fam = out_dir/"crossmapped/{ancestry}_crossmapped_plink.fam",
        inbed = out_dir/"crossmapped/{ancestry}_crossmap_input.bed",
        outbed = out_dir/"crossmapped/{ancestry}_crossmap_output.bed",
        excluded_ids = out_dir/"crossmapped/{ancestry}_excluded_ids.txt"
    params:
        in_plink = out_dir/"subset_ancestry/{ancestry}_subset",
        out = out_dir/"crossmapped/{ancestry}_crossmapped_plink",
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    container:
        config['deps']['container']
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}
        plink2 --pfile {params.in_plink} \
            --exclude {output.excluded_ids} \
            --make-bed \
            --output-chr MT \
            --out {params.out}
        awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        pgen = out_dir/"crossmapped/{ancestry}_crossmapped_plink.bed",
        psam = out_dir/"crossmapped/{ancestry}_crossmapped_plink.bim",
        pvar = out_dir/"crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    params:
        infile = out_dir/"crossmapped/{ancestry}_crossmapped_plink",
        out = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted"
    container:
        config['deps']['container']
    shell:
        """
        plink2 --bfile {params.infile} \
            --make-bed \
            --max-alleles 2 \
            --output-chr MT \
            --out {params.out}
        """


rule harmonize_hg38:
    input:
        bed = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi"
    output:
        bed = out_dir/"harmonize_hg38/{ancestry}.bed",
        bim = out_dir/"harmonize_hg38/{ancestry}.bim",
        fam = out_dir/"harmonize_hg38/{ancestry}.fam"
    resources:
        java_mem = lambda wildcards, attempt: attempt * config["imputation"]["harmonize_hg38_java_memory"],
    params:
        infile = out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted",
        out = out_dir/"harmonize_hg38/{ancestry}",
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    container:
        config['deps']['container']
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
        bed = out_dir/"harmonize_hg38/{ancestry}.bed",
        bim = out_dir/"harmonize_hg38/{ancestry}.bim",
        fam = out_dir/"harmonize_hg38/{ancestry}.fam"
    output:
        data_vcf_gz = out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index = out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    params:
        infile = out_dir/"harmonize_hg38/{ancestry}",
        out = out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38"
    container:
        config['deps']['container']
    shell:
        """
        plink2 --bfile {params.infile} \
            --recode vcf id-paste=iid \
            --chr 1-22 \
            --out {params.out}

        bgzip {params.out}.vcf
        bcftools index {output.data_vcf_gz}
        """


rule vcf_fixref_hg38:
    input:
        fasta = config['ref']['dir'] + config['ref']['fasta'],
        vcf = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz",
        index = config['ref']['dir'] + config['ref']['vcf_dir'] + "/30x-GRCh38_NoSamplesSorted.vcf.gz.tbi",
        data_vcf = out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf = out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index = out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    container:
        config['deps']['container']
    shell:
        """
        bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} | \
        bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}

        #Index
        bcftools index {output.vcf}
        """


rule filter_preimpute_vcf:
    input:
        vcf = out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf = out_dir/"filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf = out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index = out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    params:
        maf = lambda wildcards: float(maf_df["MAF"][maf_df.Ancestry == wildcards.ancestry].values),
        missing = config["imputation"]["snp_missing_pct"],
        hwe = config["imputation"]["snp_hwe"]
    container:
        config['deps']['container']
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
        vcf = out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf = temp(out_dir/"het/{ancestry}_filtered_temp.vcf"),
        inds = out_dir/"het/{ancestry}_het_failed.inds",
        het = out_dir/"het/{ancestry}_het.het",
        passed = out_dir/"het/{ancestry}_het_passed.inds",
        passed_list = out_dir/"het/{ancestry}_het_passed.txt"
    params:
        het_base = out_dir/"het/{ancestry}_het",
        script = "/opt/SNP_imputation_1000g_hg38/Imputation/scripts/filter_het.R",
        hwe = out_dir/"hwe/{ancestry}_hwe",
        out = out_dir/"het/{ancestry}_het"
    container:
        config['deps']['container']
    shell:
        """
        gunzip -c {input.vcf} \
            | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            > {output.tmp_vcf}
        vcftools --vcf {output.tmp_vcf} --het --out {params.het_base}
        Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list = out_dir/"het/{ancestry}_het_passed.txt",
        vcf = out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        index = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["het_filter_memory"]
    params:
        hwe = out_dir/"hwe/{ancestry}_hwe",
        out = out_dir/"het_filter/{ancestry}_het_filter"
    container:
        config['deps']['container']
    shell:
        """
        bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf}

        #Index the output file
        bcftools index {output.vcf}
        """


rule calculate_missingness:
    input:
        filtered_vcf = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf = temp(out_dir/"filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss = out_dir/"filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals = out_dir/"genotype_donor_annotation/{ancestry}_individuals.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * config["imputation"]["calculate_missingness_memory"]
    params:
        out = out_dir/"filter_preimpute_vcf/{ancestry}_genotypes"
    container:
        config['deps']['container']
    shell:
        """
        gunzip -c {input.filtered_vcf} | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' > {output.tmp_vcf}

        vcftools --gzvcf {output.tmp_vcf} --missing-indv --out {params.out}

        bcftools query -l {input.filtered_vcf} >> {output.individuals}
        """


rule split_by_chr:
    input:
        filtered_vcf = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index = out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf = out_dir/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        index = out_dir/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi"
    container:
        config['deps']['container']
    shell:
        """
        bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

rule eagle_prephasing:
    input:
        vcf = out_dir/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz",
        map_file = Path(config['ref']['dir'])/config['ref']['genetic_map'],
        phasing_file = Path(config['ref']['dir'])/config['ref']['phasing_dir']/"chr{chr}.bcf"
    output:
        vcf = out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    params:
        out = out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased"
    container:
        config['deps']['container']
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
        vcf = out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file = config['ref']['dir'] + config['ref']['impute_dir'] + "/chr{chr}.m3vcf.gz"
    output:
        out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    params:
        out = out_dir/"minimac_imputed/{ancestry}_chr{chr}",
        minimac4 = "/opt/bin/minimac4",
        chunk_length = config["imputation"]["chunk_length"]
    container:
        config['deps']['container']
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
        vcfs = lambda wildcards: expand(out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", chr = chromosomes, ancestry = ancestry_subsets)
    output:
        combined = out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        ind = out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    params:
        files_begin = out_dir/"minimac_imputed/{ancestry}_chr*.dose.vcf.gz"
    container:
        config['deps']['container']
    shell:
        """
        bcftools concat -Oz {params.files_begin} > {output.combined}
        bcftools index {output.combined}
        """


rule combine_vcfs_all:
    input:
        vcfs = lambda wildcards: expand(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets)
    output:
        combined = out_dir/"vcf_all_merged/imputed_hg38.vcf.gz",
        ind = out_dir/"vcf_all_merged/imputed_hg38.vcf.gz.csi"
    container:
        config['deps']['container']
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
