import pandas as pd

rule reheader_vcf:
    input:
        vcf=config['deps']['input_vcf'],
        id_map=config['deps']['id_map']
    output:
        vcf=temp(out_dir/"vcf/reheader.hg19.vcf")
    log:
        logs/"vcf_reheader.log"
    container:
        config['containers']['bcftools_biocontainer']
    shell:
        """
        bcftools reheader -s {input.id_map} {input.vcf} -o {output.vcf} > {log} 2>&1
        """

rule run_plink:
    input:
        psam=config['deps']['input_psam'],
        vcf=rules.reheader_vcf.output.vcf
    output:
        psam=temp(out_dir/"plink/plink.psam"),
        pgen=temp(out_dir/"plink/plink.pgen"),
        pvar=temp(out_dir/"plink/plink.pvar")
    log:
        logs/"plink/run_plink.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        pgen={output.pgen}
        plink_prefix=${{pgen%.pgen}}
        plink2 \
            --vcf {input.vcf} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $plink_prefix \
            --psam {input.psam} \
            --sort-vars
        mv ${{plink_prefix}}.log {log}
        """

rule calculate_missingness:
    input:
        psam=rules.run_plink.output.psam,
        pgen=rules.run_plink.output.pgen,
        pvar=rules.run_plink.output.pvar
    output:
        smiss=temp(out_dir/"indiv_missingness/sample_missingness.smiss")
    log:
        logs/"plink/calculate_missingness.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_smiss={output.smiss}
        out_prefix=${{out_smiss%.smiss}}
        plink2 --threads {threads} \
            --pfile $in_prefix \
            --missing \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """

rule indiv_missingness:
    input:
        psam=rules.run_plink.output.psam,
        pgen=rules.run_plink.output.pgen,
        pvar=rules.run_plink.output.pvar,
        smiss=rules.calculate_missingness.output.smiss
    output:
        pgen=out_dir/"indiv_missingness/indiv_missingness.pgen",
        pvar=out_dir/"indiv_missingness/indiv_missingness.pvar",
        psam=out_dir/"indiv_missingness/indiv_missingness.psam"
    params:
        mind = config["params"]["indiv_missingness_mind"]
    log:
        logs/"plink/indiv_missingness.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_bed={output.pgen}
        out_prefix=${{out_bed%.pgen}}
        plink2 --threads {threads} \
            --pfile $in_prefix \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --mind {params.mind} \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """

rule check_sex:
    input:
        pgen=rules.indiv_missingness.output.pgen,
        pvar=rules.indiv_missingness.output.pvar,
        psam=rules.indiv_missingness.output.psam
    output:
        bed=temp(out_dir/"check_sex/check_sex.bed"),
        bim=temp(out_dir/"check_sex/check_sex.bim"),
        fam=temp(out_dir/"check_sex/check_sex.fam"),
        sexcheck=temp(out_dir/"check_sex/check_sex.sexcheck"),
        sexcheck_tsv=temp(out_dir/"check_sex/check_sex.sexcheck.tsv"),
        hh=temp(out_dir/"check_sex/check_sex.hh"),
        no=temp(out_dir/"check_sex/check_sex.nosex")
    log:
        logs/"plink/check_sex.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_bed={output.bed}
        out_prefix=${{out_bed%.bed}}

        plink2 --threads {threads} \
            --pfile $in_prefix \
            --make-bed \
            --max-alleles 2 \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}

        plink --threads {threads} \
            --bfile $out_prefix \
            --check-sex \
            --biallelic-only strict \
            --out $out_prefix
        cat ${{out_prefix}}.log >> {log}
        rm ${{out_prefix}}.log

        awk '{{$1=$1}}1' OFS="\t" {output.sexcheck} > {output.sexcheck_tsv}
        touch {output.no}
        """

rule find_common_snps:
    input:
        pvar=out_dir/"indiv_missingness/indiv_missingness.pvar",
        pvar_1000g=config['refs']['1000g']['pvar']
    output:
        snps_data=temp(out_dir/"common_snps/snps_data.tsv"),
        snps_1000g=temp(out_dir/"common_snps/snps_1000g.tsv")
    log:
        logs/"find_common_snps.log"
    shell:
        """
        exec > {log} 2>&1
        # Find the intersection of SNPs from both pvar files using awk
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.pvar} {input.pvar_1000g} > {output.snps_1000g}
        awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.pvar_1000g} {input.pvar} > {output.snps_data}
        """

rule extract_common_snps:
    input:
        pgen=rules.indiv_missingness.output.pgen,
        pvar=rules.indiv_missingness.output.pvar,
        psam=rules.indiv_missingness.output.psam,
        pgen_1000g=config['refs']['1000g']['pgen'],
        pvar_1000g=config['refs']['1000g']['pvar'],
        psam_1000g=config['refs']['1000g']['psam'],
        snps_data=out_dir/"common_snps/snps_data.tsv",
        snps_1000g=out_dir/"common_snps/snps_1000g.tsv"
    output:
        pgen=temp(out_dir/"common_snps/subset_data.pgen"),
        pvar=temp(out_dir/"common_snps/subset_data.pvar"),
        psam=temp(out_dir/"common_snps/subset_data.psam"),
        pgen_1000g=temp(out_dir/"common_snps/subset_1000g.pgen"),
        pvar_1000g=temp(out_dir/"common_snps/subset_1000g.pvar"),
        psam_1000g=temp(out_dir/"common_snps/subset_1000g.psam")
    log:
        logs/"plink/extract_common_snps.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        exec > {log} 2>&1
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        echo "Input prefix: $in_prefix"
        out_pgen={output.pgen}
        out_prefix=${{out_pgen%.pgen}}
        echo "Output prefix: $out_prefix"

        cmd="plink2 --threads {threads} \
            --pfile $in_prefix \
            --extract {input.snps_data} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $out_prefix"
        echo $cmd; eval $cmd
        rm ${{out_prefix}}.log

        in_pgen_1000g={input.pgen_1000g}
        in_prefix_1000g=${{in_pgen_1000g%.pgen}}
        echo "Input prefix: $in_prefix_1000g"
        out_pgen_1000g={output.pgen_1000g}
        out_prefix_1000g=${{out_pgen_1000g%.pgen}}
        echo "Output prefix: $out_prefix_1000g"

        cmd="plink2 --threads {threads} \
            --pfile $in_prefix_1000g \
            --extract {input.snps_1000g} \
            --make-pgen \
            --out $out_prefix_1000g"
        echo $cmd; eval $cmd
        rm ${{out_prefix_1000g}}.log
        """

rule prune_1000g:
    input:
        bed_1000g=out_dir/"common_snps/subset_1000g.pgen",
        bim_1000g=out_dir/"common_snps/subset_1000g.pvar",
        fam_1000g=out_dir/"common_snps/subset_1000g.psam",
        bim=out_dir/"common_snps/subset_data.pvar",
        bed=out_dir/"common_snps/subset_data.pgen",
        fam=out_dir/"common_snps/subset_data.psam"
    output:
        prune_in_1000t=temp(out_dir/"common_snps/subset_pruned_1000g.prune.in"),
        prune_out_1000g=temp(out_dir/"common_snps/subset_pruned_1000g.prune.out"),
        prune_out=temp(out_dir/"common_snps/subset_data.prune.out"),
        bed_1000g=temp(out_dir/"common_snps/subset_pruned_1000g.pgen"),
        bim_1000g=temp(out_dir/"common_snps/subset_pruned_1000g.pvar"),
        fam_1000g=temp(out_dir/"common_snps/subset_pruned_1000g.psam"),
        bed=temp(out_dir/"common_snps/subset_pruned_data.pgen"),
        bim=temp(out_dir/"common_snps/subset_pruned_data.pvar"),
        bim_temp=temp(out_dir/"common_snps/subset_pruned_data_temp.pvar"),
        bim_old=temp(out_dir/"common_snps/subset_pruned_data_original.pvar"),
        fam=temp(out_dir/"common_snps/subset_pruned_data.psam"),
        data_1000g_key=temp(out_dir/"common_snps/subset_pruned_data_1000g_key.txt"),
        SNPs2keep=temp(out_dir/"common_snps/SNPs2keep.txt")
    params:
        ld_prune=config['params']['ld_prune']
    log:
        logs/"plink/prune_1000g.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        exec > {log} 2>&1
        in_pgen_1000g={input.bed_1000g}
        in_prefix_1000g=${{in_pgen_1000g%.pgen}}
        out_pgen_1000g={output.bed_1000g}
        out_prefix_1000g=${{out_pgen_1000g%.pgen}}  # out_dir/"common_snps/subset_pruned_1000g

        plink2 --threads {threads} \
            --pfile $in_prefix_1000g \
            --indep-pairwise {params.ld_prune} \
            --out $out_prefix_1000g
        rm ${{out_prefix_1000g}}.log

        plink2 --threads {threads} \
            --pfile $in_prefix_1000g \
            --extract {output.prune_out_1000g} \
            --make-pgen \
            --out $out_prefix_1000g
        rm ${{out_prefix_1000g}}.log

        if [[ $(grep "##" {input.bim} | wc -l) > 0 ]]
        then
            grep "##" {input.bim} > {output.data_1000g_key}
        fi
        awk -F"\\t" 'BEGIN{{OFS=FS = "\\t"}} \
            NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} \
            {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.bim_1000g} {input.bim} \
            | grep -v "##" >> {output.data_1000g_key}
        grep -v "##" {output.data_1000g_key} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' \
            > {output.prune_out}

        in_pgen={input.bed}
        in_prefix=${{in_pgen%.pgen}}
        out_pgen={output.bed}
        out_prefix=${{out_pgen%.pgen}}
        plink2 --threads {threads} \
            --pfile $in_prefix \
            --extract {output.prune_out} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $out_prefix
        rm ${{out_prefix}}.log

        cp {output.bim} {output.bim_old}
        grep -v "#" {output.bim_old} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' \
            > {output.SNPs2keep}
        grep "#CHROM" {output.data_1000g_key} > {output.bim}
        grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        grep "##" {output.bim_1000g} > {output.bim}
        cat {output.bim_temp} >> {output.bim}
        """

rule final_pruning: ### put in contingency for duplicated snps - remove from both 1000G and your dataset
    input:
        bed=out_dir/"common_snps/subset_pruned_data.pgen",
        bim=out_dir/"common_snps/subset_pruned_data.pvar",
        fam=out_dir/"common_snps/subset_pruned_data.psam"
    output:
        bed=temp(out_dir/"common_snps/final_subset_pruned_data.pgen"),
        bim=temp(out_dir/"common_snps/final_subset_pruned_data.pvar"),
        fam=temp(out_dir/"common_snps/final_subset_pruned_data.psam")
    log:
        logs/"plink/final_pruning.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.bed}
        in_prefix=${{in_pgen%.pgen}}
        out_pgen={output.bed}
        out_prefix=${{out_pgen%.pgen}}
        plink2 --rm-dup 'force-first' \
            --threads {threads} \
            --pfile $in_prefix \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """

# KING-robust kinship is ancestry-agnostic, so this runs on the pooled pruned
# data to catch cross-ancestry duplicates and sample swaps. Two plink2 passes
# are needed because --king-cutoff prunes samples before --make-king-table
# writes the .kin0, so combining them would silently truncate the report.
# Upstream pfiles all set psam-cols=fid,..., so .king.cutoff.out.id is
# guaranteed two-column FID/IID and the tail-append in rule separate_indivs
# is safe.
rule relatedness_check:
    input:
        bed=rules.final_pruning.output.bed,
        bim=rules.final_pruning.output.bim,
        fam=rules.final_pruning.output.fam
    output:
        kin0=out_dir/"relatedness/relatedness_check.kin0",
        remove_id=temp(out_dir/"relatedness/relatedness_check.king.cutoff.out.id"),
        in_id=temp(out_dir/"relatedness/relatedness_check.king.cutoff.in.id")
    params:
        king_cutoff=config["params"]["king_cutoff"],
        king_table_cutoff=config["params"]["king_table_cutoff"]
    log:
        logs/"plink/relatedness_check.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.bed}
        in_prefix=${{in_pgen%.pgen}}
        out_kin0={output.kin0}
        out_prefix=${{out_kin0%.kin0}}

        plink2 --threads {threads} \
            --pfile $in_prefix \
            --make-king-table \
            --king-table-filter {params.king_table_cutoff} \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}

        plink2 --threads {threads} \
            --pfile $in_prefix \
            --king-cutoff {params.king_cutoff} \
            --out $out_prefix
        cat ${{out_prefix}}.log >> {log}
        rm ${{out_prefix}}.log
        """

### use PCA from plink for PCA and projection
rule pca_1000g:
    input:
        pgen_1000g=out_dir/"common_snps/subset_pruned_1000g.pgen",
        pvar_1000g=out_dir/"common_snps/subset_pruned_1000g.pvar",
        psam_1000g=out_dir/"common_snps/subset_pruned_1000g.psam",
        bed=out_dir/"common_snps/subset_pruned_data.pgen"
    output:
        out=temp(out_dir/"pca_projection/subset_pruned_1000g_pcs.acount"),
        eig_all=temp(out_dir/"pca_projection/subset_pruned_1000g_pcs.eigenvec.allele"),
        eig_vec=temp(out_dir/"pca_projection/subset_pruned_1000g_pcs.eigenvec"),
        eig=temp(out_dir/"pca_projection/subset_pruned_1000g_pcs.eigenval")
    log:
        logs/"plink/pca_1000g.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen_1000g}
        in_prefix=${{in_pgen%.pgen}}
        out_eig={output.eig}
        out_prefix=${{out_eig%.eigenval}}
        plink2 --threads {threads} \
            --pfile $in_prefix \
            --freq counts \
            --pca allele-wts \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """

### use plink pca results to plot with R ###
rule pca_project:
    input:
        pgen=out_dir/"common_snps/final_subset_pruned_data.pgen",
        pvar=out_dir/"common_snps/final_subset_pruned_data.pvar",
        psam=out_dir/"common_snps/final_subset_pruned_data.psam",
        pgen_1000g=out_dir/"common_snps/subset_pruned_1000g.pgen",
        pvar_1000g=out_dir/"common_snps/subset_pruned_1000g.pvar",
        psam_1000g=out_dir/"common_snps/subset_pruned_1000g.psam",
        frq=out_dir/"pca_projection/subset_pruned_1000g_pcs.acount",
        scores=out_dir/"pca_projection/subset_pruned_1000g_pcs.eigenvec.allele"
    output:
        projected_scores=temp(out_dir/"pca_projection/final_subset_pruned_data_pcs.sscore"),
        projected_1000g_scores=temp(out_dir/"pca_projection/subset_pruned_1000g_pcs_projected.sscore")
    log:
        logs/"plink/pca_project.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_scores={output.projected_scores}
        out_prefix=${{out_scores%.sscore}}
        plink2 --threads {threads} \
            --pfile $in_prefix \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}

        in_pgen_1000g={input.pgen_1000g}
        in_prefix_1000g=${{in_pgen_1000g%.pgen}}
        out_scores_1000g={output.projected_1000g_scores}
        out_prefix_1000g=${{out_scores_1000g%.sscore}}
        plink2 --threads {threads} \
            --pfile $in_prefix_1000g \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation variance-standardize \
            --score-col-nums 6-15 \
            --out $out_prefix_1000g
        cat ${{out_prefix_1000g}}.log >> {log}
        rm ${{out_prefix_1000g}}.log
       """

rule pca_projection_assign:
    input:
        projected_scores = out_dir/"pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = out_dir/"pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = out_dir/"common_snps/subset_1000g.psam",
        psam = out_dir/"plink/plink.psam",
        sexcheck = out_dir/"check_sex/check_sex.sexcheck.tsv"
    output:
        sexcheck=temp(out_dir/"pca_sex_checks/check_sex_update_remove.tsv"),
        anc_check=temp(out_dir/"pca_sex_checks/ancestry_update_remove.tsv"),
        plot=out_dir/"pca_sex_checks/Ancestry_PCAs.png"
    container:
        config['containers']['r']
    # conda:
    #     "../envs/r.yaml"
    log:
        logs/"pca_projection_assign.log"
    script:
        "../scripts/PCA_Projection_Plotting.R"

rule separate_indivs:
    input:
        sexcheck=out_dir/"pca_sex_checks/check_sex_update_remove.tsv",
        anc_check=out_dir/"pca_sex_checks/ancestry_update_remove.tsv",
        king_remove=rules.relatedness_check.output.remove_id
    output:
        update_sex=temp(out_dir/"separate_indivs/sex_update_indivs.tsv"),
        remove_indiv=temp(out_dir/"separate_indivs/remove_indivs.tsv"),
        remove_indiv_temp=temp(out_dir/"separate_indivs/remove_indivs_temp.tsv")
    log:
        logs/"plink/separate_indivs.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        exec > {log} 2>&1
        grep "UPDATE" {input.sexcheck} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2,$4)}}' \
            | sed 's/SNPSEX/SEX/g'\
            > {output.update_sex}
        grep "REMOVE" {input.sexcheck} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}'\
            > {output.remove_indiv_temp}
        grep "REMOVE" {input.anc_check} \
            | awk 'BEGIN{{FS=OFS="\t"}}{{print($1,$2)}}' \
            >> {output.remove_indiv_temp}
        tail -n +2 {input.king_remove} \
            >> {output.remove_indiv_temp}
        sort -u {output.remove_indiv_temp} > {output.remove_indiv}
        """

# This rule updates the psam file to include ancestry
rule update_psam:
    input:
        ancestry=out_dir/"pca_sex_checks/ancestry_update_remove.tsv",
        sex=out_dir/"pca_sex_checks/check_sex_update_remove.tsv",
        psam=out_dir/"indiv_missingness/indiv_missingness.psam"
    output:
        anc_updated_psam=temp(out_dir/"pca_sex_checks/updated_psam.psam")
    log:
        logs/"update_psam.log"
    run:
        # Read input files
        ancestry_check = pd.read_csv(input.ancestry, sep="\t")
        psam_df_local = pd.read_csv(input.psam, sep="\t")
        # Filter the ancestry check for individuals to be updated
        ancestry_updates = ancestry_check[ancestry_check["UPDATE/REMOVE/KEEP"] == "UPDATE"]
        # Merge the PSAM DataFrame with the ancestry updates based on 'IID'
        updated_psam = psam_df_local.merge(
            ancestry_updates[["IID", "PCA_Assignment"]],
            on="IID",
            how="left"
        )
        # Update 'Provided_Ancestry' where ancestry updates are available
        updated_psam["Provided_Ancestry"] = updated_psam["PCA_Assignment"].combine_first(updated_psam["Provided_Ancestry"])
        # Drop the 'PCA_Assignment' column used for the merge
        updated_psam.drop(columns=["PCA_Assignment"], inplace=True)
        # Write updated PSAM to output file
        updated_psam.to_csv(output.anc_updated_psam, sep="\t", na_rep="NA", index=False)

rule prepare_update:
    input:
        pgen=out_dir/"indiv_missingness/indiv_missingness.pgen",
        pvar=out_dir/"indiv_missingness/indiv_missingness.pvar",
        anc_updated_psam=out_dir/"pca_sex_checks/updated_psam.psam"
    output:
        pgen=temp(out_dir/"update_prep/update_prep.pgen"),
        pvar=temp(out_dir/"update_prep/update_prep.pvar"),
        psam=temp(out_dir/"update_prep/update_prep.psam")
    log:
        logs/"prepare_update.log"
    shell:
        """
        cp {input.pgen} {output.pgen}
        cp {input.pvar} {output.pvar}
        cp {input.anc_updated_psam} {output.psam}
        """

rule update_sex_ancestry:
    input:
        pgen=out_dir/"update_prep/update_prep.pgen",
        pvar=out_dir/"update_prep/update_prep.pvar",
        psam=out_dir/"update_prep/update_prep.psam",
        update_sex=out_dir/"separate_indivs/sex_update_indivs.tsv",
        remove_indiv=out_dir/"separate_indivs/remove_indivs.tsv",
    output:
        pgen=final/"post_qc.pgen",
        pvar=final/"post_qc.pvar",
        psam=final/"post_qc.psam"
    log:
        logs/"update_sex_ancestry.log"
    container:
        config['containers']['plink']
    # conda:
    #     "../envs/plink.yaml"
    shell:
        """
        in_pgen={input.pgen}
        in_prefix=${{in_pgen%.pgen}}
        out_pgen={output.pgen}
        out_prefix=${{out_pgen%.pgen}}

        plink2 --threads {threads} \
            --pfile $in_prefix \
            --update-sex {input.update_sex} \
            --remove {input.remove_indiv} \
            --make-pgen 'psam-cols='fid,parents,sex,phenos \
            --out $out_prefix
        mv ${{out_prefix}}.log {log}
        """

rule export_sample_metadata:
    input:
        psam=final/"post_qc.psam"
    output:
        csv=metadata_csv
    log:
        logs/"export_sample_metadata.log"
    run:
        psam_df = pd.read_csv(input.psam, sep="\t", dtype=str)

        # Select relevant columns
        df = psam_df[["IID", "SEX", "Provided_Ancestry"]].copy()

        # Map sex codes to labels
        sex_map = {"1": "Male", "2": "Female"}
        df["SEX"] = df["SEX"].map(sex_map).fillna("Unknown")

        # HANCESTRO sample-level inferred ancestry
        ancestry_map = {
            "EUR": "European ancestry",
            "EAS": "East Asian ancestry",
            "AFR": "African ancestry",
            "SAS": "South Asian ancestry",
            "AMR": "Latin American or Admixed American ancestry",
        }
        ancestry_id_map = {
            "EUR": "HANCESTRO:0005",
            "EAS": "HANCESTRO:0009",
            "AFR": "HANCESTRO:0010",
            "SAS": "HANCESTRO:0006",
            "AMR": "HANCESTRO:0014",
        }

        # HANCESTRO 1KGP reference superpopulation (provenance)
        ref_superpop_id_map = {
            "EUR": "HANCESTRO:2003",
            "EAS": "HANCESTRO:2002",
            "AFR": "HANCESTRO:2000",
            "SAS": "HANCESTRO:2004",
            "AMR": "HANCESTRO:2001",
        }

        raw_anc = df["Provided_Ancestry"]
        df["ancestry"] = raw_anc.map(ancestry_map)
        df["ancestry_id"] = raw_anc.map(ancestry_id_map)
        df["reference_superpopulation"] = raw_anc.where(raw_anc.isin(ref_superpop_id_map), other=raw_anc)
        df["reference_superpopulation_id"] = raw_anc.map(ref_superpop_id_map)

        # Rename and select final columns
        df = df.rename(columns={"IID": "gussid", "SEX": "sex"})
        df = df[["gussid", "sex", "ancestry", "ancestry_id", "reference_superpopulation", "reference_superpopulation_id"]]

        df.to_csv(output.csv, index=False)

