rule subset_psam_by_ancestry:
    input:
        psam = out_dir/"update_sex_ancestry/update_sex.psam"
    output:
        keep=temp(out_dir/"subset_ancestry/{ancestry}_individuals.psam")
    shell:
        """
        grep {wildcards.ancestry} {input.psam} > {output.keep}
        """

rule subset_plink_by_ancestry:
    input:
        keep=out_dir/"subset_ancestry/{ancestry}_individuals.psam",
        pgen=out_dir/"update_sex_ancestry/update_sex.pgen",
        psam=out_dir/"update_sex_ancestry/update_sex.psam",
        pvar=out_dir/"update_sex_ancestry/update_sex.pvar"
    output:
        pgen = out_dir/"subset_ancestry/{ancestry}_subset.pgen",
        psam = out_dir/"subset_ancestry/{ancestry}_subset.psam",
        pvar = out_dir/"subset_ancestry/{ancestry}_subset.pvar"
    log:
        logs/"subset_plink_by_ancestry_{ancestry}.log"
    container:
        config['deps']['container']
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_pgen={output.pgen}
        out_prefix=${{out_pgen%.pgen}}

        plink2 --threads {threads} \
            --pfile $in_prefix \
            --keep {input.keep}  \
            --max-alleles 2 \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $out_prefix \
            > {log} 2>&1
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
        chain_file = "/opt/GRCh37_to_GRCh38.chain"
    container:
        config['deps']['container']
    shell:
        """
        awk 'BEGIN{{FS=OFS="\t"}}{{print $1,$2,$2+1,$3,$4,$5}}' {input.pvar} > {output.inbed}
        CrossMap.py bed {params.chain_file} {output.inbed} {output.outbed}
        awk '{{print $4}}' {output.outbed}.unmap > {output.excluded_ids}

        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_bed={output.bed}
        out_prefix=${{out_bed%.bed}}

        plink2 --pfile $in_prefix \
            --exclude {output.excluded_ids} \
            --make-bed \
            --output-chr MT \
            --out $out_prefix

        awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        bed=out_dir/"crossmapped/{ancestry}_crossmapped_plink.bed",
        bim=out_dir/"crossmapped/{ancestry}_crossmapped_plink.bim",
        fam=out_dir/"crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam"
    log:
        logs/"sort_bed_{ancestry}.log"
    container:
        config['deps']['container']
    shell:
        """
        in_bed={input.bed}
        in_prefix=${{in_bed%.bed}}
        out_bed={output.bed}
        out_prefix=${{out_bed%.bed}}

        plink2 --bfile $in_prefix \
            --make-bed \
            --max-alleles 2 \
            --output-chr MT \
            --out $out_prefix \
            > {log} 2>&1
        """


rule harmonize_hg38:
    input:
        bed=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf=config['refs']['vcf'],
        index=config['refs']['vcf'] + ".tbi"
    output:
        bed=out_dir/"harmonize_hg38/{ancestry}.bed",
        bim=out_dir/"harmonize_hg38/{ancestry}.bim",
        fam=out_dir/"harmonize_hg38/{ancestry}.fam"
    params:
        jar = "/opt/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    log:
        logs/"harmonize_hg38_{ancestry}.log"
    container:
        config['deps']['container']
    shell:
        """
        in_bed={input.bed}
        in_prefix=${{in_bed%.bed}}
        out_bed={output.bed}
        out_prefix=${{out_bed%.bed}}

        java -Xmx{resources.java_mem}g -jar {params.jar}\
            --input $in_prefix \
            --inputType PLINK_BED \
            --ref {input.vcf} \
            --refType VCF \
            --update-id \
            --output $out_prefix
            > {log} 2>&1
        """

rule plink_to_vcf:
    input:
        bed=out_dir/"harmonize_hg38/{ancestry}.bed",
        bim=out_dir/"harmonize_hg38/{ancestry}.bim",
        fam=out_dir/"harmonize_hg38/{ancestry}.fam"
    output:
        vcf=out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz",
        index=out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi"
    log:
        logs/"plink_to_vcf_{ancestry}.log"
    container:
        config['deps']['container']
    shell:
        """
        in_bed={input.bed}
        in_prefix=${{in_bed%.bed}}
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.vcf.gz}}

        plink2 --bfile $in_prefix \
            --recode vcf id-paste=iid \
            --chr 1-22 \
            --out $out_prefix \
            > {log} 2>&1

        bgzip ${{out_prefix}}.vcf
        bcftools index {output.vcf}
        """

rule vcf_fixref_hg38:
    input:
        fasta=config['refs']['hg38_fa'],
        vcf=config['refs']['vcf'],
        index=config['refs']['vcf'] + ".tbi",
        data_vcf=out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf=out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz",
        index=out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi"
    container:
        config['deps']['container']
    shell:
        """
        bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} \
            | bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

rule filter_preimpute_vcf:
    input:
        vcf=out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf=out_dir/"filter_preimpute_vcf/{ancestry}_tagged.vcf.gz",
        filtered_vcf=out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
        filtered_index=out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi"
    params:
        maf=config['params']['maf'],
        missing=config["params"]["snp_missing_pct"],
        hwe=config["params"]["snp_hwe"]
    container:
        config['deps']['container']
    shell:
        """
        #Add tags
        export BCFTOOLS_PLUGINS=/opt/bcftools-1.10.2/plugins
        bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} \
            | bcftools filter -e 'REF="N" | REF="I" | REF="D"' \
            | bcftools filter -e "ALT='.'" \
            | bcftools norm -d all \
            | bcftools norm -m+any \
            | bcftools view -m2 -M2 -Oz -o {output.filtered_vcf}

        #Index the output file
        bcftools index {output.filtered_vcf}
        """

rule het:
    input:
        vcf=out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf=temp(out_dir/"het/{ancestry}_filtered_temp.vcf"),
        inds=out_dir/"het/{ancestry}_het_failed.inds",
        het=out_dir/"het/{ancestry}_het.het",
        passed=out_dir/"het/{ancestry}_het_passed.inds",
        passed_list=out_dir/"het/{ancestry}_het_passed.txt"
    params:
        script="/opt/SNP_imputation_1000g_hg38/Imputation/scripts/filter_het.R"
    container:
        config['deps']['container']
    shell:
        """
        het={output.het}
        het_base=${{het%.het}}
        gunzip -c {input.vcf} \
            | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            > {output.tmp_vcf}
        vcftools --vcf {output.tmp_vcf} --het --out $het_base
        Rscript {params.script} {output.het} {output.inds} {output.passed} {output.passed_list}
        """

rule het_filter:
    input:
        passed_list=out_dir/"het/{ancestry}_het_passed.txt",
        vcf=out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"
    output:
        vcf=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        index=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
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
        filtered_vcf=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf=temp(out_dir/"filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss=out_dir/"filter_preimpute_vcf/{ancestry}_genotypes.imiss",
        individuals=out_dir/"genotype_donor_annotation/{ancestry}_individuals.tsv"
    container:
        config['deps']['container']
    shell:
        """
        gunzip -c {input.filtered_vcf} \
            | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            > {output.tmp_vcf}

        out_miss={output.miss}
        out_prefix=${{out_miss%.imiss}}
        vcftools --gzvcf {output.tmp_vcf} --missing-indv --out $out_prefix

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
        map_file = config['refs']['genetic_map'],
        phasing_file = config['refs']['phasing'] + "chr{chr}.bcf"
    output:
        vcf = out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz"
    container:
        config['deps']['container']
    shell:
        """
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.vcf.gz}}
        eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix=$out_prefix \
            --numThreads={threads}
        """

rule minimac_imputation:
    input:
        vcf=out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file=config['refs']['impute'] + "/chr{chr}.m3vcf.gz"
    output:
        vcf=out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"
    params:
        minimac4 = "/opt/bin/minimac4",
        chunk_length = config["params"]["chunk_length"]
    container:
        config['deps']['container']
    shell:
        """
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.dose.vcf.gz}}
        {params.minimac4} --refHaps {input.impute_file} \
            --haps {input.vcf} \
            --prefix $out_prefix \
            --format GT,DS,GP \
            --noPhoneHome \
            --cpus {threads} \
            --ChunkLengthMb {params.chunk_length}
        """

rule combine_vcfs_ancestry:
    input:
        vcfs = lambda wildcards: expand(out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", ancestry=wildcards.ancestry, chr=chromosomes)
    output:
        combined = out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz",
        ind = out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi"
    log:
        logs/"combine_vcfs_{ancestry}.log"
    container:
        config['deps']['container']
    shell:
        """
        bcftools concat -Oz {input.vcfs} > {output.combined} 2> {log}
        bcftools index {output.combined} >> {log} 2>&1
        """

rule combine_vcfs_all:
    input:
        vcfs = expand(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets)
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

rule restore_vcf_header:
    input:
        vcf=out_dir/"vcf_all_merged/imputed_hg38.vcf.gz",
        id_map=config['deps']['id_map']
    output:
        final=out_dir/"imputed_hg38.vcf"
    log:
        logs/"restore_vcf_header.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """
        awk -F'\t' '{{print $2 "\t" $1}}' {input.id_map} \
            | bcftools reheader -s - {input.vcf} -o {output.final} \
            2> {log}
        bcftools sort {output.final} 2>> {log}
        """