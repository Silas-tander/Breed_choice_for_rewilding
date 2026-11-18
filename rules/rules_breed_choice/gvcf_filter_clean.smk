# Site and genotype filters of the raw joint-called gVCFs (per autosome), with bcftools stats per step.

# Curated RefSeq accession -> chr1-31 lookup used by the rename rules
rule stage_chr_names_breed_choice:
    input:
        "data/chr"
    output:
        directory("results/breed_choice/chr")
    log:
        "logs/breed_choice_setup/stage_chr_names_breed_choice.log"
    threads: 1
    resources:
        mem_mb = 500,
        runtime = 10
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        mkdir -p {output}
        cp {input}/* {output}/
        """

rule bcf_stats_raw:
    input:
        vcf = "data/gvcf/{autosome}.gvcf.gz"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_raw.stats"
    threads: 2
    resources:
        mem_mb = 16000,
        runtime = 680
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule bcf_stats_var:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.bcf.filt.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}.var.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule plot_gvcf_stats:
    input:
        path = "data/gvcf/{autosome}.gvcf.gz",
    params:
        outdire = "results/breed_choice/plots/{autosome}/"
    output:
        DP_plot = "results/breed_choice/plots/{autosome}/Histogram_of_Depth_Values_Stratified_by_SNP_Status.png"
    threads: 1
    resources:
        mem_mb = 24000,
        runtime = 480
    conda:
        "../../envs/slurm_v1.yaml"
    shell:
        "python scripts/python/gwf_plot_diagnostics.py {input.path} {params.outdire}"

rule bcftools_harsh_filter:
    input:
        vcf = "data/gvcf/{autosome}.gvcf.gz",
    params:
        qual = 20,  # not used: the expression has no QUAL filter
        dp = 15*146,  # INFO/DP window: 15-30x per sample, 146 samples
        dp_max = 30*146
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        bcftools view -i '({params.dp_max} > INFO/DP) & (INFO/DP > {params.dp}) & (INFO/AC > 0) & (INFO/ExcessHet >= 1E-6)' {input.vcf} -Ov -o {output.vcf_filt}
        '''

rule bcf_stats_harsh:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_bcf_harsh.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule remove_genmap:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf",
        bed_mask = "data/ref_genome/genmap_mask.merge.bed",
        harsh_stats = "results/breed_choice/genomics_stats/filt/{autosome}_bcf_harsh.stats"  # ordering only: stats read the VCF before it is deleted
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.genmap.filt.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 720
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        vcftools --gzvcf {input.vcf} \
        --exclude-bed {input.bed_mask} \
        --recode --stdout > {output.vcf_filt}
        rm {input.vcf}
        '''

rule remove_indels:
    input:
        vcf="results/breed_choice/gvcf_filt/{autosome}.genmap.filt.gvcf"
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.final.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        bcftools view {input.vcf} --types snps -M 2 -m 2 > {output.vcf_filt}
        rm {input.vcf}
        '''

rule bcf_stats:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.final.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_final.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule remove_repeats:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.final.gvcf",
        bed_mask = "data/ref_genome/repeats_Equcab3.bed",
        final_stats = "results/breed_choice/genomics_stats/filt/{autosome}_final.stats",  # ordering only: readers finish before the VCF is deleted
        gtmiss_wreps = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.reps.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 720
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        vcftools --vcf {input.vcf} \
        --exclude-bed {input.bed_mask} \
        --recode --stdout > {output.vcf_filt}
        rm {input.vcf}
        '''

rule bcf_stats_reps:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.reps.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}.reps.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

# Genotypes -> missing: DP < 8 or > 90, allele balance < 0.05 (hets), GQ < 8; sites are kept
rule bcf_setGT:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.reps.gvcf"
    params:
        AD_threshold = 0.05
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        # Compress the input VCF with bgzip if not already compressed
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        
        # Apply the filtering steps using bcftools
        bcftools view {input.vcf}.gz | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP<8' | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP>90' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:1] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:0] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/GQ<8' | \
        bcftools view -o {output.vcf_filt}
        
        # Clean up intermediate files if desired
        rm {input.vcf}.gz {input.vcf}.gz.tbi
        '''


# Same with repeats kept; input of heterozygosity and SnpEff genetic load
rule bcf_setGT_with_repeats:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.final.gvcf"
    params:
        AD_threshold = 0.05
    output:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        # Compress the input VCF with bgzip if not already compressed
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz

        # Apply the filtering steps using bcftools
        bcftools view {input.vcf}.gz | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP<8' | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP>90' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:1] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:0] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/GQ<8' | \
        bcftools view -o {output.vcf_filt}

        # Clean up intermediate files if desired
        rm {input.vcf}.gz {input.vcf}.gz.tbi
        '''

rule bcf_stats_GTmiss:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf",
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_GTmiss.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule bcf_maf_filt:
    input:
        vcf_filt = "results/breed_choice/gvcf_filt/{autosome}.reps.gvcf"
    output:
        vcf_maf = "results/breed_choice/gvcf_filt/{autosome}.final.maf.gvcf.gz"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        '''
        bcftools view {input.vcf_filt} -q 0.05:minor -Oz -o {output.vcf_maf}
        bcftools index {output.vcf_maf}
        '''

# LD_pruning_vcf / LD_pruning_plink are alternatives; LD_pruning_bcf is the one used downstream
rule LD_pruning_vcf:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/breed_choice/gvcf_filt/{autosome}_GTmiss.LDprune.gvcf.recode.vcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "vcftools --vcf {input.vcf} --thin 300000 --recode --out {output.vcf}"

rule LD_pruning_plink:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_plink.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/plink.yml"
    shell:
        '''
        plink --indep-pairwise 50 100 0.2 --out {input.vcf} --allow-extra-chr --chr-set 31 no-xy no-mt && \
        plink --extract {input.vcf}.prune.in --recode vcf --out {output.vcf} --allow-extra-chr --chr-set 31 no-xy no-mt
        '''

rule LD_pruning_bcf:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools +prune -l 0.2 -w 50000 {input.vcf} -o {output.vcf}"

rule gather_split_vcf_LD:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf",
            autosome=config["autosomes_only"]
        )
    output:
        "results/breed_choice/gathered_vcf/all_LD_pruned_bcf.gvcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -o {output}
        """

rule bcf_stats_maf:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.final.maf.gvcf.gz"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_maf.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule bcf_stats_GTmiss_wReps:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_GTmiss_wReps.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"
