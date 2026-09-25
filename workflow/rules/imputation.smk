rule subset_psam_by_ancestry:
    input:
        psam=post_qc_plink/"post_qc.psam"
    output:
        keep=temp(out_dir/"subset_ancestry/{ancestry}_individuals.psam")
    shell:
        """
        grep {wildcards.ancestry} {input.psam} > {output.keep}
        """

rule subset_plink_by_ancestry:
    input:
        keep=out_dir/"subset_ancestry/{ancestry}_individuals.psam",
        pgen=post_qc_plink/"post_qc.pgen",
        psam=post_qc_plink/"post_qc.psam",
        pvar=post_qc_plink/"post_qc.pvar"
    output:
        pgen=temp(out_dir/"subset_ancestry/{ancestry}_subset.pgen"),
        psam=temp(out_dir/"subset_ancestry/{ancestry}_subset.psam"),
        pvar=temp(out_dir/"subset_ancestry/{ancestry}_subset.pvar")
    log:
        logs/"subset_plink_by_ancestry/subset_{ancestry}.log"
    conda:
        "../envs/crossmap.yaml"
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
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
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
        bed=temp(out_dir/"crossmapped/{ancestry}_crossmapped_plink.bed"),
        bim=temp(out_dir/"crossmapped/{ancestry}_crossmapped_plink.bim"),
        fam=temp(out_dir/"crossmapped/{ancestry}_crossmapped_plink.fam"),
        inbed=temp(out_dir/"crossmapped/{ancestry}_crossmap_input.bed"),
        outbed=temp(out_dir/"crossmapped/{ancestry}_crossmap_output.bed"),
        excluded_ids=temp(out_dir/"crossmapped/{ancestry}_excluded_ids.txt"),
        unmap=temp(out_dir/"crossmapped/{ancestry}_crossmap_output.bed.unmap")
    params:
        chain_file = "resources/liftover/GRCh37_to_GRCh38.chain.gz"
    conda:
        "../envs/crossmap.yaml"
    log:
        logs/"crossmap/crossmap_{ancestry}.log"
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
        mv ${{out_prefix}}.log {log}

        awk -F'\t' 'BEGIN {{OFS=FS}} {{print $1,$4,0,$2,$6,$5}}' {output.outbed} > {output.bim}
        """

rule sort_bed:
    input:
        bed=out_dir/"crossmapped/{ancestry}_crossmapped_plink.bed",
        bim=out_dir/"crossmapped/{ancestry}_crossmapped_plink.bim",
        fam=out_dir/"crossmapped/{ancestry}_crossmapped_plink.fam"
    output:
        bed=temp(out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed"),
        bim=temp(out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim"),
        fam=temp(out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam")
    log:
        logs/"sort_bed/sort_bed_{ancestry}.log"
    conda:
        "../envs/plink.yaml"
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
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """


rule harmonize_hg38:
    input:
        bed=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bed",
        bim=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.bim",
        fam=out_dir/"crossmapped_sorted/{ancestry}_crossmapped_sorted.fam",
        vcf=config['refs']['vcf'],
        index=config['refs']['vcf'] + ".tbi"
    output:
        bed=temp(out_dir/"harmonize_hg38/{ancestry}.bed"),
        bim=temp(out_dir/"harmonize_hg38/{ancestry}.bim"),
        fam=temp(out_dir/"harmonize_hg38/{ancestry}.fam"),
        updates=temp(out_dir/"harmonize_hg38/{ancestry}_idUpdates.txt")
    params:
        jar = "resources/tools/GenotypeHarmonizer-1.4.23/GenotypeHarmonizer.jar"
    log:
        harmonizer=logs/"harmonize_hg38/harmonize_hg38_{ancestry}.log",
        snp_log=logs/"harmonize_hg38/snpLog_{ancestry}.log"
    conda:
        "../envs/java.yaml"
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
        mv $(dirname $out_bed)/{wildcards.ancestry}.log {log.harmonizer}
        mv $(dirname $out_bed)/{wildcards.ancestry}_snpLog.log {log.snp_log}
        """

rule plink_to_vcf:
    input:
        bed=out_dir/"harmonize_hg38/{ancestry}.bed",
        bim=out_dir/"harmonize_hg38/{ancestry}.bim",
        fam=out_dir/"harmonize_hg38/{ancestry}.fam"
    output:
        vcf=temp(out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"),
        index=temp(out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz.csi")
    log:
        logs/"plink_to_vcf_{ancestry}.log"
    conda:
        "../envs/plink-bcftools.yaml"
    shell:
        """
        in_bed={input.bed}
        in_prefix=${{in_bed%.bed}}
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.vcf.gz}}

        plink2 --bfile $in_prefix \
            --recode vcf id-paste=iid \
            --chr 1-22 \
            --out $out_prefix

        mv ${{out_prefix}}.log {log}

        bgzip ${{out_prefix}}.vcf
        bcftools index {output.vcf}
        """

rule vcf_fixref_hg38:
    input:
        fasta=config['refs']['hg38_int_fa'],
        vcf=config['refs']['vcf'],
        index=config['refs']['vcf'] + ".tbi",
        data_vcf=out_dir/"harmonize_hg38/{ancestry}_harmonised_hg38.vcf.gz"
    output:
        vcf=temp(out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"),
        index=temp(out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz.csi")
    log:
        logs/"vcf_fixref_hg38/fixref_{ancestry}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools +fixref {input.data_vcf} -- -f {input.fasta} -i {input.vcf} \
            | bcftools norm --check-ref x -f {input.fasta} -Oz -o {output.vcf} \
            2> {log}
        bcftools index {output.vcf} 2>> {log}
        """

rule filter_preimpute_vcf:
    input:
        vcf=out_dir/"vcf_fixref_hg38/{ancestry}_fixref_hg38.vcf.gz"
    output:
        tagged_vcf=temp(out_dir/"filter_preimpute_vcf/{ancestry}_tagged.vcf.gz"),
        filtered_vcf=temp(out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz"),
        filtered_index=temp(out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz.csi")
    params:
        maf=config['params']['pre_maf'],
        missing=config["params"]["snp_missing_pct"],
        hwe=config["params"]["snp_hwe"]
    log:
        logs/"filter_preimpute_vcf/filter_{ancestry}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        #Add tags
        export BCFTOOLS_PLUGINS="$CONDA_PREFIX/libexec/bcftools"
        bcftools +fill-tags {input.vcf} -Oz -o {output.tagged_vcf} 2> {log}

        #Filter rare and non-HWE variants and those with abnormal alleles and duplicates
        bcftools filter -i 'INFO/HWE > {params.hwe} & F_MISSING < {params.missing} & MAF[0] > {params.maf}' {output.tagged_vcf} \
            | bcftools filter -e 'REF="N" | REF="I" | REF="D"' \
            | bcftools filter -e "ALT='.'" \
            | bcftools norm -d all \
            | bcftools norm -m+any \
            | bcftools view -m2 -M2 -Oz -o {output.filtered_vcf} \
            2>> {log}

        #Index the output file
        bcftools index {output.filtered_vcf} 2>> {log}
        """

rule het:
    input:
        vcf=out_dir/"filter_preimpute_vcf/{ancestry}_filtered.vcf.gz",
    output:
        tmp_vcf=temp(out_dir/"het/{ancestry}_filtered_temp.vcf"),
        inds=temp(out_dir/"het/{ancestry}_het_failed.inds"),
        het=temp(out_dir/"het/{ancestry}_het.het"),
        passed=temp(out_dir/"het/{ancestry}_het_passed.inds"),
        passed_list=temp(out_dir/"het/{ancestry}_het_passed.txt")
    params:
        script=workflow.source_path("../scripts/filter_het.R")
    conda:
        "../envs/het.yaml"
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
        vcf=temp(out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz"),
        index=temp(out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi")
    log:
        logs/"het_filter/het_filter_{ancestry}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -S {input.passed_list} {input.vcf} -Oz -o {output.vcf} 2> {log}

        #Index the output file
        bcftools index {output.vcf} 2>> {log}
        """

rule calculate_missingness:
    input:
        filtered_vcf=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        tmp_vcf=temp(out_dir/"filter_preimpute_vcf/{ancestry}_het_filtered.vcf"),
        miss=temp(out_dir/"filter_preimpute_vcf/{ancestry}_genotypes.imiss"),
        individuals=temp(out_dir/"genotype_donor_annotation/{ancestry}_individuals.tsv")
    log:
        logs/"calculate_missingness/missingness_{ancestry}.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        gunzip -c {input.filtered_vcf} \
            | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
            > {output.tmp_vcf}

        out_miss={output.miss}
        out_prefix=${{out_miss%.imiss}}
        vcftools --gzvcf {output.tmp_vcf} --missing-indv --out $out_prefix 2> {log}

        bcftools query -l {input.filtered_vcf} >> {output.individuals} 2>> {log}
        """

rule split_by_chr:
    input:
        filtered_vcf=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz",
        filtered_index=out_dir/"het_filter/{ancestry}_het_filtered.vcf.gz.csi"
    output:
        vcf=temp(out_dir/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz"),
        index=temp(out_dir/"split_by_chr/{ancestry}_chr_{chr}.vcf.gz.csi")
    log:
        logs/"split_by_chr/split_{ancestry}_chr{chr}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools view -r {wildcards.chr} {input.filtered_vcf} -Oz -o {output.vcf} 2> {log}
        bcftools index {output.vcf} 2>> {log} 
        """

rule eagle_prephasing:
    input:
        vcf=rules.split_by_chr.output.vcf,
        index=rules.split_by_chr.output.index,
        map_file = config['refs']['genetic_map'],
        phasing_file = config['refs']['phasing'] + "chr{chr}.bcf"
    output:
        vcf=temp(out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz")
    log:
        logs/"eagle/eagle_prephasing_{ancestry}_chr{chr}.log"
    conda:
        "../envs/eagle.yaml"
    shell:
        """
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.vcf.gz}}
        eagle --vcfTarget={input.vcf} \
            --vcfRef={input.phasing_file} \
            --geneticMapFile={input.map_file} \
            --chrom={wildcards.chr} \
            --outPrefix=$out_prefix \
            --numThreads={threads} \
            > {log} 2>&1
        """

rule minimac_imputation:
    input:
        vcf=out_dir/"eagle_prephasing/{ancestry}_chr{chr}_phased.vcf.gz",
        impute_file=config['refs']['impute'] + "/chr{chr}.m3vcf.gz"
    output:
        vcf=temp(out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz"),
        info=temp(out_dir/"minimac_imputed/{ancestry}_chr{chr}.info")
    params:
        minimac4 = "resources/tools/minimac4",
        chunk_length = config["params"]["chunk_length"]
    log:
        logs/"minimac/minimac_{ancestry}_chr{chr}.log"
    conda:
        "../envs/minimac4.yaml"
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
            --ChunkLengthMb {params.chunk_length} \
            > {log} 2>&1
        """

rule combine_vcfs_ancestry:
    input:
        vcfs=lambda wildcards: expand(out_dir/"minimac_imputed/{ancestry}_chr{chr}.dose.vcf.gz", ancestry=wildcards.ancestry, chr=chromosomes)
    output:
        combined=temp(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz"),
        ind=temp(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi")
    log:
        logs/"combine_vcfs_{ancestry}.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        bcftools concat -Oz {input.vcfs} > {output.combined} 2> {log}
        bcftools index {output.combined} >> {log} 2>&1
        """

rule combine_vcfs_all:
    input:
        vcfs=expand(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz", ancestry = ancestry_subsets),
        indices=expand(out_dir/"vcf_merged_by_ancestries/{ancestry}_imputed_hg38.vcf.gz.csi", ancestry = ancestry_subsets)
    output:
        combined=temp(out_dir/"vcf_all_merged/imputed_hg38.vcf.gz"),
        ind=temp(out_dir/"vcf_all_merged/imputed_hg38.vcf.gz.csi")
    log:
        logs/"combine_vcfs_all.log"
    conda:
        "../envs/bcftools.yaml"
    shell:
        """
        if [[ $(ls -l {input.vcfs} | wc -l) > 1 ]]
        then
            bcftools merge -Oz {input.vcfs} > {output.combined} 2> {log}
        else
            cp {input.vcfs} {output.combined}
        fi
        bcftools index {output.combined} >> {log} 2>&1
        """

rule restore_vcf_header:
    input:
        vcf=out_dir/"vcf_all_merged/imputed_hg38.vcf.gz",
        id_map=config['deps']['id_map']
    output:
        reheadered=temp(out_dir/"restore_header/imputed_hg38.vcf.gz"),
        idx=temp(out_dir/"restore_header/imputed_hg38.vcf.gz.csi")
    log:
        logs/"restore_vcf_header.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """
        awk -F'\t' '{{print $2 "\t" $1}}' {input.id_map} \
            | bcftools reheader -s - {input.vcf} -o {output.reheadered} \
            2> {log}
        bcftools index {output.reheadered} 2>> {log}
        """

# NOTE: The input VCF (pre-QC) may contain samples that were removed during
# upstream QC steps (e.g. individual missingness filtering via --mind in the
# plink_qc stage). We must subset the XY chromosomes to only post-QC samples
# so that reinsert_XY can merge them with the imputed autosomes without a
# sample mismatch.
rule get_preimpute_XY:
    input:
        preimpute=config['deps']['input_vcf'],
        map=config['refs']['chr2int'],
        psam=post_qc_plink/"post_qc.psam"
    output:
        input_idx=temp(config['deps']['input_vcf']+'.csi'),
        pre_XY=temp(out_dir/"preimpute_XY.vcf.gz"),
        pre_XY_idx=temp(out_dir/"preimpute_XY.vcf.gz.csi"),
        samples=temp(out_dir/"post_qc_samples.txt")
    log:
        logs/"preimpute_XY.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """
        awk 'NR>1 {{print $2}}' {input.psam} > {output.samples}
        bcftools index {input.preimpute}
        bcftools view -r chrX,chrY -S {output.samples} {input.preimpute} \
            | bcftools annotate --rename-chrs {input.map} \
            | bcftools sort -Oz -o {output.pre_XY} \
            2> {log}
        bcftools index {output.pre_XY} 2>> {log}
        """

rule reinsert_XY:
    input:
        xy=rules.get_preimpute_XY.output.pre_XY,
        xy_idx=rules.get_preimpute_XY.output.pre_XY_idx,
        auto=rules.restore_vcf_header.output.reheadered,
        auto_idx=rules.restore_vcf_header.output.idx
    output:
        order=temp(out_dir/"sample_order.txt"),
        reordered_auto=temp(out_dir/"reordered_auto.vcf.gz"),
        vcf=temp(out_dir/"merged/imputed_hg38.vcf.gz"),
        index=temp(out_dir/"merged/imputed_hg38.vcf.gz.csi")
    log:
        logs/"reinsert_XY.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """
        bcftools query -l {input.xy} > {output.order}
        bcftools view -S {output.order} -o {output.reordered_auto} {input.auto}
        bcftools concat {output.reordered_auto} {input.xy} -Oz -o {output.vcf} 2> {log}
        bcftools index {output.vcf}  2>> {log}
        """

rule filter_maf_r2:
    input:
        vcf=rules.reinsert_XY.output.vcf,
    output:
        vcf=temp(out_dir/"filtered/imputed_filtered_maf_r2.hg38.vcf.gz"),
        idx=temp(out_dir/"filtered/imputed_filtered_maf_r2.hg38.vcf.gz.csi")
    params:
        maf=config['params']['post_maf']
    log:
        logs/"filter_maf_r2.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """ 
            bcftools filter -i '(IMPUTED=1 && MAF >= {params.maf} && R2 > 0.8) || (IMPUTED=0)' \
                -Oz -o {output.vcf} \
                {input.vcf} \
                2> {log}
            bcftools index {output.vcf} 2>> {log}
        """

rule filter_exons_indels:
    input:
        vcf=rules.filter_maf_r2.output.vcf,
        bed=config['refs']['bed']
    output:
        vcf=temp(out_dir/"filter_exons_indels/imputed_filtered_maf_r2.hg38.recode.vcf.gz")
    log:
        logs/"filter_exons.indels.log"
    conda:
        "../envs/vcftools.yaml"
    shell:
        """
        out_vcf={output.vcf}
        out_prefix=${{out_vcf%.recode.vcf.gz}}
        vcftools --gzvcf {input.vcf} \
            --max-alleles 2 \
            --remove-indels \
            --bed {input.bed} \
            --recode \
            --recode-INFO-all \
            --out $out_prefix \
            2> {log}
        bgzip ${{out_prefix}}.recode.vcf 2>> {log}
        """

rule rename_chromosomes:
    input:
        vcf=rules.filter_exons_indels.output.vcf,
        map=config['refs']['int2chr'],
        fai=config['refs']['hg38_chr_fai']
    output:
        tmp_vcf=out_dir/"temp.vcf",
        vcf=out_dir/"imputed_filtered.hg38.vcf.gz"
    log:
        logs/"rename_chromosomes.log"
    container:
        "docker://quay.io/biocontainers/bcftools:1.21--h8b25389_0"
    shell:
        """
        bcftools annotate --rename-chrs {input.map} {input.vcf} \
            | bcftools view \
            | grep -v "##contig=<ID=" \
            > {output.tmp_vcf}
        bcftools reheader --fai {input.fai} {output.tmp_vcf} \
            | bcftools sort \
            | sed 's/; Date=.*//g' \
            | bgzip -c \
            > {output.vcf}
        bcftools index {output.vcf}
        """
