
# rule PLINK_ped_map:
#     input:
#         vcf = "results/variants/bcftools_call/{sample}_chrom.bcf",
#     params:
#         MAF = 0.05
#     output:
#         # directory("results/variants/plink_files/{sample}/")
#         path = "results/variants/plink_files/",
#         mapa = "results/variants/plink_files/{sample}.map",
#         ped = "results/variants/plink_files/{sample}.ped",
#         nosex = "results/variants/plink_files/{sample}.nosex"
#     threads: 2
#     resources:
#         mem_mb = 1000,
#         runtime = 120
#     conda:
#         "../envs/plink.yml"
#     shell:
#         "plink --bfile {input.vcf} --recode --out {output.path}"

# rule gather_split_gvcf:
#     input:
#         gvcfs = expand("results/variants/genomics_filt_vcf/{interval}.filt.vcf",
#                         interval = config["intervals"])
#     output:
#         "results/variants/gathered_filt/{sample}.all.filt.g.vcf"
#     threads: 4
#     resources:
#         mem_mb = 20000,
#         runtime = 120
#     conda:
#         "../envs/gatk.yaml"
#     script:
#         "../scripts/python/gather_gvcf.py"

# rule change_chr_name:
#     input:
#         vcf = "results/variants/genomics_filt_vcf/{interval}.filt.vcf"
#     shell:
#         "bcftools annotate --rename-chrs chr_name_conv.txt original.vcf.gz | bgzip > rename.vcf.gz"

rule get_chr2:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.final.gvcf",
    output:
        chrs = "results/gwf/chr/{autosome}.chr.txt"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools query -f '%CHROM\n' {input.vcf}|uniq > {output.chrs}"

rule rename_chr2:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.gvcf",
    params:
        chrs = "results/gwf/chr/{autosome}.chr.txt",
        names = "results/gwf/chr/chr.names.txt"
    output:
        vcf_re = "results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.rename.gvcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re}
        # rm {input.vcf}
        '''

rule PLINK_ped_map_vcf2:
    input:
        vcf = "results/gwf/gvcf_filt/{autosome}.GTmiss_wReps.rename.gvcf"
    params:
        MAF = 0.05,
        path = "results/gwf/plink_files.GTmiss_wReps/{autosome}"
    output:
        # directory("results/variants/plink_files/{sample}/")
        mapa = "results/gwf/plink_files.GTmiss_wReps/{autosome}.map",
        ped = "results/gwf/plink_files.GTmiss_wReps/{autosome}.ped",
        nosex = "results/gwf/plink_files.GTmiss_wReps/{autosome}.nosex"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/plink.yml"
    shell:
        "plink --vcf {input.vcf} --recode --allow-extra-chr --chr-set 31 no-xy no-mt --out {params.path}"

rule gather_split_vcf_shell:
    input:
        vcf_ints = expand("results/gwf/gvcf_filt/{autosome}.final.gvcf",
                    autosome=config["autosomes"])
    output:
        "results/gwf/gathered_vcf/vcfs/all_final.vcf"
    threads: 4
    resources:
        mem_mb=12000,
        runtime=720 // 2
    conda:
        "../envs/gatk.yaml"
    shell:
        '''
        gatk GatherVcfs \
            --INPUT {input.vcf_ints} \
            --OUTPUT {output}
        '''

rule rezipping:
    input:
        vcf = "results/snpeff/gathered_vcf/{sample}.ann.vcf",
    output:
        bgzipped = "results/snpeff/gathered_vcf/{sample}.ann.vcf.gz"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../envs/bcftools.yml"
    shell:
        '''
        zcat {input.vcf} | bgzip -c > {output.bgzipped} && rm {input.vcf}
        '''

# merge for different samples, concat for same samples concat should be the answer here
rule gather_split_vcf_all:
    input:
        vcf_ints=expand(
            "results/snpeff/gathered_vcf/{sample}.ann.vcf.gz",
            sample=config["samples"]
        )
    output:
        merged_vcf="results/snpeff/gathered_vcf/all/all_samples.ann.vcf.gz",
        merged_index="results/snpeff/gathered_vcf/all/all_samples.ann.vcf.gz.csi"
    threads: 4
    resources:
        mem_mb=4*8000,
        runtime=720 // 2
    conda:
        "../envs/bcftools.yml"
    shell:
        """
        # Ensure all input VCFs are indexed before concatenation
        for vcf in {input.vcf_ints}; do
            if [ ! -f "$vcf.csi" ]; then
                bcftools index "$vcf"
            fi
        done

        # Concatenate indexed VCFs
        bcftools concat -Oz -o {output.merged_vcf} {input.vcf_ints}

        # Index the merged VCF
        bcftools index {output.merged_vcf}
        """

rule get_chr_all2:
    input:
        vcf = "results/gwf/gathered_vcf/all_final.vcf"
    output:
        chrs = "results/gwf/chr/all.chr.txt"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools query -f '%CHROM\n' {input.vcf} | uniq > {output.chrs}"

rule rename_chr_all2:
    input:
        vcf = "results/gwf/gathered_vcf/vcfs/all_final.vcf"
    params:
        # chrs = "results/gwf/chr/{autosome}.chr.txt",
        names = "results/gwf/chr/chr.names.txt"
    output:
        vcf_re = "results/gwf/gathered_vcf/vcfs/all_final.rename.vcf"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../envs/bcftools.yml"
    shell:
        "bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re}"

rule PLINK_ped_map_all2:
    input:
        vcf = "results/gwf/gathered_vcf/vcfs/all_final.rename.vcf"
    params:
        MAF = 0.05,
        path = "results/gwf/gathered_vcf/plink/all_final"
    output:
        # directory("results/variants/plink_files/{sample}/")
        mapa = "results/gwf/gathered_vcf/plink/all_final.map",
        ped = "results/gwf/gathered_vcf/plink/all_final.ped",
        nosex = "results/gwf/gathered_vcf/plink/all_final.nosex"
    threads: 2
    resources:
        mem_mb = 1000,
        runtime = 120
    conda:
        "../envs/plink.yml"
    shell:
        "plink --vcf {input.vcf} --recode --allow-extra-chr --chr-set 31 no-xy no-mt --out {params.path}"

rule PLINK_vcf_bed:
    input:
        vcf = "results/gwf/gathered_vcf/all_LD_pruned_bcf.gvcf"
    params:
        path = "results/gwf/gathered_vcf/plink/LD_pruned",
        bim = "results/gwf/gathered_vcf/plink/LD_pruned.bim",
    output:
        bed = "results/gwf/gathered_vcf/plink/LD_pruned.bed",
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../envs/plink.yml"
    shell:
        """
        plink --vcf {input.vcf} --allow-extra-chr --const-fid --geno 0.999 --make-bed --out {params.path}
        """

#### Do this to bim file to not have chromosomal naming problems
# awk '{$1="0";print $0}' {params.bim} > {params.bim}.tmp
# mv {params.bim}.tmp {params.bim}
