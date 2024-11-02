#!/usr/bin/env python
shell.executable('bash')

rule indiv_missingness:
    input:
        pgen = plink/"plink.pgen",
        pvar = plink/"plink.pvar",
        psam = plink/"plink.psam",
    output:
        bed = imputation/"indiv_missingness/indiv_missingness.pgen",
        bim = imputation/"indiv_missingness/indiv_missingness.pvar",
        fam = imputation/"indiv_missingness/indiv_missingness.psam",
    params:
       infile = plink/"plink",
       out = imputation/"indiv_missingness/indiv_missingness",
       mind = config["plink_gender_ancestry_QC"]["indiv_missingness_mind"]
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        echo {params.infile}
        plink2 --threads {threads} \
            --pfile {params.infile} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --mind {params.mind} \
            --out {params.out}
        """

## Updated on 6 June by Drew to include --biallelic-only strict
## Updated on 7 June by Drew Neavin to remove from output: hh = imputation/"check_sex/check_sex.hh",
rule check_sex:
    input:
        bed = imputation/"indiv_missingness/indiv_missingness.pgen",
        bim = imputation/"indiv_missingness/indiv_missingness.pvar",
        fam = imputation/"indiv_missingness/indiv_missingness.psam",
    output:
        bed = imputation/"check_sex/check_sex.bed",
        bim = imputation/"check_sex/check_sex.bim",
        fam = imputation/"check_sex/check_sex.fam",
        log = imputation/"check_sex/check_sex.log",
        nosex = imputation/"check_sex/check_sex.nosex",
        sexcheck = imputation/"check_sex/check_sex.sexcheck",
        sexcheck_tab = imputation/"check_sex/check_sex.sexcheck.tsv"
    params:
        infile = imputation/"indiv_missingness/indiv_missingness",
        out = imputation/"check_sex/check_sex"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile} --make-bed --max-alleles 2 --out {params.out}
        plink --threads {threads} --bfile {params.out} --check-sex --biallelic-only strict --out {params.out}
        touch {output.nosex}
        sed 's/^ \+//g' {output.sexcheck} | sed 's/ \+/\t/g' > {output.sexcheck_tab}
        """

### Pull just common SNPs between two groups ###
rule common_snps:
    input:
        bed = imputation/"indiv_missingness/indiv_missingness.pgen",
        bim = imputation/"indiv_missingness/indiv_missingness.pvar",
        fam = imputation/"indiv_missingness/indiv_missingness.psam",
    output:
        snps_data = imputation/"common_snps/snps_data.tsv",
        snps_1000g = imputation/"common_snps/snps_1000g.tsv",
        bed = imputation/"common_snps/subset_data.pgen",
        bim = imputation/"common_snps/subset_data.pvar",
        fam = imputation/"common_snps/subset_data.psam",
        bed_1000g = imputation/"common_snps/subset_1000g.pgen",
        bim_1000g = imputation/"common_snps/subset_1000g.pvar",
        fam_1000g = imputation/"common_snps/subset_1000g.psam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["common_snps_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * config["plink_gender_ancestry_QC"]["common_snps_memory"]
    params:
        bim_1000 = "/opt/1000G/all_phase3_filtered.pvar",
        infile = imputation/"indiv_missingness/indiv_missingness",
        infile_1000g = "/opt/1000G/all_phase3_filtered",
        out = imputation/"common_snps/subset_data",
        out_1000g = imputation/"common_snps/subset_1000g"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim} {params.bim_1000} > {output.snps_1000g}
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {params.bim_1000} {input.bim} > {output.snps_data}
        plink2 --threads {threads} --pfile {params.infile} --extract {output.snps_data} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.snps_1000g} --make-pgen --out {params.out_1000g}
        """

### Prune with --indep,
rule prune_1000g:
    input:
        bed_1000g = imputation/"common_snps/subset_1000g.pgen",
        bim_1000g = imputation/"common_snps/subset_1000g.pvar",
        fam_1000g = imputation/"common_snps/subset_1000g.psam",
        bim = imputation/"common_snps/subset_data.pvar",
        bed = imputation/"common_snps/subset_data.pgen",
        fam = imputation/"common_snps/subset_data.psam",
    output:
        prune_out_1000g = imputation/"common_snps/subset_pruned_1000g.prune.out",
        prune_out = imputation/"common_snps/subset_data.prune.out",
        bed_1000g = imputation/"common_snps/subset_pruned_1000g.pgen",
        bim_1000g = imputation/"common_snps/subset_pruned_1000g.pvar",
        fam_1000g = imputation/"common_snps/subset_pruned_1000g.psam",
        bed = imputation/"common_snps/subset_pruned_data.pgen",
        bim = imputation/"common_snps/subset_pruned_data.pvar",
        bim_temp = imputation/"common_snps/subset_pruned_data_temp.pvar",
        bim_old = imputation/"common_snps/subset_pruned_data_original.pvar",
        fam = imputation/"common_snps/subset_pruned_data.psam",
        data_1000g_key = imputation/"common_snps/subset_pruned_data_1000g_key.txt",
        SNPs2keep = imputation/"common_snps/SNPs2keep.txt"
    params:
        out_1000g = imputation/"common_snps/subset_pruned_1000g",
        infile_1000g = imputation/"common_snps/subset_1000g",
        infile = imputation/"common_snps/subset_data",
        out = imputation/"common_snps/subset_pruned_data"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile_1000g} \
            --indep-pairwise 50 5 0.5 \
            --out {params.out_1000g}
        plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.prune_out_1000g} --make-pgen --out {params.out_1000g}
        if [[ $(grep "##" {input.bim} | wc -l) > 0 ]]
        then
            grep "##" {input.bim} > {output.data_1000g_key}
        fi
        awk -F"\\t" 'BEGIN{{OFS=FS = "\\t"}} NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.bim_1000g} {input.bim} | grep -v "##" >> {output.data_1000g_key}
        grep -v "##" {output.data_1000g_key} | awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' > {output.prune_out}
        plink2 --threads {threads} --pfile {params.infile} --extract {output.prune_out} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        cp {output.bim} {output.bim_old}
        grep -v "#" {output.bim_old} | awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        grep "#CHROM" {output.data_1000g_key} > {output.bim}
        grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        grep "##" {output.bim_1000g} > {output.bim}
        cat {output.bim_temp} >> {output.bim}
        """
        
rule final_pruning: ### put in contingency for duplicated snps - remove from both 1000G and your dataset
    input:
        bed = imputation/"common_snps/subset_pruned_data.pgen",
        bim = imputation/"common_snps/subset_pruned_data.pvar",
        fam = imputation/"common_snps/subset_pruned_data.psam",
    output:
        bed = imputation/"common_snps/final_subset_pruned_data.pgen",
        bim = imputation/"common_snps/final_subset_pruned_data.pvar",
        fam = imputation/"common_snps/final_subset_pruned_data.psam",
    params:
        infile = imputation/"common_snps/subset_pruned_data",
        out = imputation/"common_snps/final_subset_pruned_data"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --rm-dup 'force-first' -threads {threads} --pfile {params.infile} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """


### use PCA from plink for PCA and projection
rule pca_1000g:
    input:
        bed_1000g = imputation/"common_snps/subset_pruned_1000g.pgen",
        bim_1000g = imputation/"common_snps/subset_pruned_1000g.pvar",
        fam_1000g = imputation/"common_snps/subset_pruned_1000g.psam",
        bed = imputation/"common_snps/subset_pruned_data.pgen" 
    output:
        out = imputation/"pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = imputation/"pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = imputation/"pca_projection/subset_pruned_1000g_pcs.eigenvec",
        eig = imputation/"pca_projection/subset_pruned_1000g_pcs.eigenval",
    params:
        infile = imputation/"common_snps/subset_pruned_1000g",
        out = imputation/"pca_projection/subset_pruned_1000g_pcs"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        plink2 --threads {threads} --pfile {params.infile} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


### use plink pca results to plot with R ###
rule pca_project:
    input:
        bed = imputation/"common_snps/final_subset_pruned_data.pgen",
        bim = imputation/"common_snps/final_subset_pruned_data.pvar",
        fam = imputation/"common_snps/final_subset_pruned_data.psam",
        frq = imputation/"pca_projection/subset_pruned_1000g_pcs.acount",
        scores = imputation/"pca_projection/subset_pruned_1000g_pcs.eigenvec.allele"
    output:
        projected_scores = imputation/"pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = imputation/"pca_projection/subset_pruned_1000g_pcs_projected.sscore"
    params:
        infile = imputation/"common_snps/final_subset_pruned_data",
        infile_1000g = imputation/"common_snps/subset_pruned_1000g",
        out = imputation/"pca_projection/final_subset_pruned_data_pcs",
        out_1000g = imputation/"pca_projection/subset_pruned_1000g_pcs_projected"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        export OMP_NUM_THREADS={threads}
        plink2 --threads {threads} --pfile {params.infile} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        plink2 --threads {threads} --pfile {params.infile_1000g} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """

rule pca_projection_assign:
    input:
        projected_scores = imputation/"pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = imputation/"pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = imputation/"common_snps/subset_1000g.psam",
        psam = plink/"plink.psam",
        sexcheck = imputation/"check_sex/check_sex.sexcheck.tsv",
    output:
        sexcheck = imputation/"pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = imputation/"pca_sex_checks/ancestry_update_remove.tsv",
    params:
        outdir = imputation/"pca_sex_checks/",
        script = "/opt/SNP_imputation_1000g_hg38/Imputation/scripts/PCA_Projection_Plotting.R"
    container:
        config['containers']['SNP_imputation']
    log:
        logs/"imputation/pca_sex_checks/variables.tsv"
    shell:
        """
        echo {params.outdir} > {log}
        echo {input.projected_scores} >> {log}
        echo {input.projected_1000g_scores} >> {log}
        echo {input.fam_1000g} >> {log}
        echo {input.psam} >> {log}
        echo {input.sexcheck} >> {log}
        Rscript {params.script} {log}
        # Mark ancestry and sex files for update
        sed -i '1! s/$/UPDATE/' {output.anc_check}
        sed -i '1! s/$/UPDATE/' {output.sexcheck}
        """

rule separate_indivs:
    input:
        sexcheck = imputation/"pca_sex_checks/check_sex_update_remove.tsv",
        anc_check = imputation/"pca_sex_checks/ancestry_update_remove.tsv"
    output:
        update_sex = imputation/"separate_indivs/sex_update_indivs.tsv",
        remove_indiv = imputation/"separate_indivs/remove_indivs.tsv",
        remove_indiv_temp = imputation/"separate_indivs/remove_indivs_temp.tsv",
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        grep "UPDATE" {input.sexcheck} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$4)}}' | sed 's/SNPSEX/SEX/g' > {output.update_sex}
        grep "REMOVE" {input.sexcheck} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}'> {output.remove_indiv_temp}
        grep "REMOVE" {input.anc_check} | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' >> {output.remove_indiv_temp}
        sort -u {output.remove_indiv_temp} > {output.remove_indiv}
        """

rule update_sex_ancestry:
    input:
        bim = imputation/"indiv_missingness/indiv_missingness.pgen",
        psam = plink/"plink.psam",
        update_sex = imputation/"separate_indivs/sex_update_indivs.tsv",
        remove_indiv = imputation/"separate_indivs/remove_indivs.tsv",
    output:
        bed = imputation/"update_sex_ancestry/update_sex.pgen",
        bim = imputation/"update_sex_ancestry/update_sex.pvar",
        psam = imputation/"update_sex_ancestry/update_sex.psam",
    params:
        anc_updated_psam = imputation/"pca_sex_checks/updated_psam.psam",
        infile = imputation/"indiv_missingness/indiv_missingness",
        psam_temp = imputation/"update_sex_ancestry/temp/indiv_missingness.psam_temp",
        tdir = imputation/"update_sex_ancestry/temp/",
        out = imputation/"update_sex_ancestry/update_sex"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        mkdir -p {params.tdir}
        cp {params.infile}* {params.tdir}
        cp {params.anc_updated_psam} {params.tdir}/indiv_missingness.psam 
        plink2 --threads {threads} --pfile {params.tdir}/indiv_missingness --update-sex {input.update_sex} --remove {input.remove_indiv} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """


rule subset_ancestry:
    input:
        psam = imputation/"update_sex_ancestry/update_sex.psam"
    output:
        keep = imputation/"subset_ancestry/{ancestry}_individuals.psam",
        pgen = imputation/"subset_ancestry/{ancestry}_subset.pgen",
        psam = imputation/"subset_ancestry/{ancestry}_subset.psam",
        pvar = imputation/"subset_ancestry/{ancestry}_subset.pvar"
    params:
        infile = imputation/"update_sex_ancestry/update_sex",
        out = imputation/"subset_ancestry/{ancestry}_subset"
    container:
        config['containers']['SNP_imputation']
    shell:
        """
        grep {wildcards.ancestry} {input.psam} > {output.keep}
        plink2 --threads {threads} --pfile {params.infile} --keep {output.keep}  --max-alleles 2 --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        """
