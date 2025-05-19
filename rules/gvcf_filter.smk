rule bcftools_harsh_filter:
    input:
        vcf = "results/gwf/gvcf/{autosome}/{autosome}.gvcf.gz",
    params:
        qual = 20,
        dp = 15*165,
        dp_max = 30*165
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        bcftools view -i '({params.dp_max} > INFO/DP) & (INFO/DP > {params.dp}) & (INFO/AC > 0) & (INFO/ExcessHet >= 1E-6)' {input.vcf} -Ov -o {output.vcf_filt}
        '''

rule bcf_stats_harsh:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf"
    output:
        vcf = "results/gwf/genomics_stats/filt/{autosome}_bcf_harsh.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule remove_genmap:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.bcf.harsh.filt.gvcf",
        bed_mask = "data/ref_genome/genmap_mask.merge.bed"
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.genmap.filt.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 720
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        vcftools --gzvcf {input.vcf} \
        --exclude-bed {input.bed_mask} \
        --recode --stdout > {output.vcf_filt}
        rm {input.vcf}
        '''

rule remove_indels:
    input:
        vcf="results/gwf/gvcf_filt/{autosome}.genmap.filt.gvcf"
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.final.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        bcftools view {input.vcf} --types snps -M 2 -m 2 > {output.vcf_filt}
        rm {input.vcf}
        '''

rule remove_repeats:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.final.gvcf",
        bed_mask = "data/ref_genome/repeats_Equcab3.bed"
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.reps.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 720
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        vcftools --vcf {input.vcf} \
        --exclude-bed {input.bed_mask} \
        --recode --stdout > {output.vcf_filt}
        rm {input.vcf}
        '''

rule bcf_stats_reps:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.reps.gvcf"
    output:
        vcf = "results/gwf/genomics_stats/filt/{autosome}.reps.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

# filtering for minor allele frequencies per individual where allele depth of minor allele is below 5%
rule bcf_setGT:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.reps.gvcf"
    params:
        AD_threshold = 0.05
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        # Compress the input VCF with bgzip if not already compressed
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        
        # Apply the filtering steps using bcftools
        bcftools view {input.vcf}.gz | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP<8' | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP>39' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:1] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:0] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/GQ<8' | \
        bcftools view -o {output.vcf_filt}
        
        # Clean up intermediate files if desired
        rm {input.vcf}.gz {input.vcf}.gz.tbi
        '''

# filtering for minor allele frequencies per individual where allele depth of minor allele is below 5%
rule bcf_setGT_with_repeats:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.final.gvcf"
    params:
        AD_threshold = 0.05
    output:
        vcf_filt = "results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        # Compress the input VCF with bgzip if not already compressed
        bgzip -c {input.vcf} > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        
        # Apply the filtering steps using bcftools
        bcftools view {input.vcf}.gz | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP<8' | \
        bcftools +setGT -- -t q -n . -i 'FMT/DP>39' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:1] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/AD[:0] / (FMT/AD[:1] + FMT/AD[:0]) < {params.AD_threshold} & FMT/GT="0/1"' | \
        bcftools +setGT -- -t q -n . -i 'FMT/GQ<8' | \
        bcftools view -o {output.vcf_filt}
        
        # Clean up intermediate files if desired
        rm {input.vcf}.gz {input.vcf}.gz.tbi
        '''

rule bcf_stats_GTmiss:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/gwf/genomics_stats/filt/{autosome}_GTmiss.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule LD_pruning_vcf:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/gwf/gvcf_filt/{autosome}_GTmiss.LDprune.gvcf.recode.vcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "vcftools --vcf {input.vcf} --thin 300000 --recode --out {output.vcf}"

rule LD_pruning_plink:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/gwf/gvcf_LD/{autosome}_GTmiss.LDprune_plink.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/plink.yml"
    shell:
        '''
        plink --indep-pairwise 50 100 0.2 --out {input.vcf} --allow-extra-chr --chr-set 31 no-xy no-mt && \
        plink --extract {input.vcf}.prune.in --recode vcf --out {output.vcf} --allow-extra-chr --chr-set 31 no-xy no-mt
        '''

rule LD_pruning_bcf:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss.gvcf"
    output:
        vcf = "results/gwf/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools +prune -l 0.2 -w 50000 {input.vcf} -o {output.vcf}"

# merge for different samples, concat for same samples concat should be the answer here
rule gather_split_vcf_LD:
    input:
        vcf_ints=lambda wildcards: expand(
            "results/gwf/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf",
            autosome=config["autosomes"]
        )
    output:
        "results/gwf/gathered_vcf/all_LD_pruned_bcf.gvcf"
    threads: 2
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -o {output}
        """

rule bcf_stats_maf:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.final.maf.gvcf.gz"
    output:
        vcf = "results/gwf/genomics_stats/filt/{autosome}_maf.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"

rule TsTv:
    input:
        vcf="results/variants/genomics_novo_vcf/{interval}.filt.SNPs.vcf"
    params:
        path = "results/variants/genomics_filt_stats/TsTv_count/{interval}"
    output:
        tstv = "results/variants/genomics_filt_stats/TsTv_count/{interval}.TsTv.count"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 30
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        vcftools --vcf {input.vcf} \
        --TsTv-by-count \
        --out {params.path}
        '''

rule TsTv_summary:
    input:
        vcf="results/variants/genomics_novo_vcf/{interval}.filt.SNPs.vcf"
    params:
        path = "results/variants/genomics_filt_stats/TsTv_summary/{interval}"
    output:
        tstv = "results/variants/genomics_filt_stats/TsTv_summary/{interval}.TsTv.summary"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        vcftools --vcf {input.vcf} \
        --TsTv-summary \
        --out {params.path}
        '''