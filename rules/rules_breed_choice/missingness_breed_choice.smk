# Per-sample and per-site missingness of the full cohort before (repeat-masked) and after the
# GTmiss genotype filters, plus a 50k-SNP genotype matrix for co-missingness.

_MISS_BC_DIR       = "results/breed_choice/genomics_stats/missingness"
_MISS_BC_PLINK_DIR = "results/breed_choice/gathered_vcf/plink_missing_breed_choice"
_MISS_BC_NAMES     = "results/breed_choice/chr/chr.names.txt"
_MISS_BC_CHR_SET   = 31

rule gather_split_vcf_reps_breed_choice:
    input:
        vcf_ints = expand(
            "results/breed_choice/gvcf_filt/{autosome}.reps.gvcf",
            autosome=config["autosomes_only"]
        )
    output:
        "results/breed_choice/gathered_vcf/all_reps.gvcf"
    log:
        "logs/breed_choice_missingness/gather_split_vcf_reps_breed_choice.log"
    threads: 2
    resources:
        mem_mb  = 12000,
        runtime = 3*720
    conda:
        "../../envs/bcftools.yml"
    shell:
        """
        bcftools concat {input.vcf_ints} -O z -o {output} 2> {log}
        """

# Raw calls use RefSeq accessions as CHROM; rename to 1-31
rule rename_chr_GTmiss_breed_choice:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.gvcf",
        chr_dir = "results/breed_choice/chr"  # ordering only: needs stage_chr_names_breed_choice
    params:
        names = _MISS_BC_NAMES
    output:
        vcf_re = "results/breed_choice/gathered_vcf/all_GTmiss.rename.gvcf"
    log:
        "logs/breed_choice_missingness/rename_chr_GTmiss_breed_choice.log"
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re} 2> {log}
        bcftools index --csi {output.vcf_re} 2>> {log}
        """

rule rename_chr_reps_breed_choice:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_reps.gvcf",
        chr_dir = "results/breed_choice/chr"  # ordering only: needs stage_chr_names_breed_choice
    params:
        names = _MISS_BC_NAMES
    output:
        vcf_re = "results/breed_choice/gathered_vcf/all_reps.rename.gvcf"
    log:
        "logs/breed_choice_missingness/rename_chr_reps_breed_choice.log"
    threads: 1
    resources:
        mem_mb  = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re} 2> {log}
        bcftools index --csi {output.vcf_re} 2>> {log}
        """

# Full cohort, before any relatedness or breed exclusion
rule plink_make_bed_GTmiss_breed_choice:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.rename.gvcf"
    output:
        bed = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bed",
        bim = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bim",
        fam = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.fam"
    params:
        base = f"{_MISS_BC_PLINK_DIR}/all_GTmiss"
    log:
        "logs/breed_choice_missingness/plink_make_bed_GTmiss_breed_choice.log"
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 180
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {_MISS_BC_PLINK_DIR}
        plink --vcf {input.vcf} \
              --allow-extra-chr --double-id \
              --chr-set {_MISS_BC_CHR_SET} no-xy no-mt \
              --set-missing-var-ids @:# \
              --make-bed \
              --threads {threads} \
              --out {params.base} \
              >> {log} 2>&1
        """

rule plink_make_bed_reps_breed_choice:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_reps.rename.gvcf"
    output:
        bed = f"{_MISS_BC_PLINK_DIR}/all_reps.bed",
        bim = f"{_MISS_BC_PLINK_DIR}/all_reps.bim",
        fam = f"{_MISS_BC_PLINK_DIR}/all_reps.fam"
    params:
        base = f"{_MISS_BC_PLINK_DIR}/all_reps"
    log:
        "logs/breed_choice_missingness/plink_make_bed_reps_breed_choice.log"
    threads: 2
    resources:
        mem_mb  = 8000,
        runtime = 180
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {_MISS_BC_PLINK_DIR}
        plink --vcf {input.vcf} \
              --allow-extra-chr --double-id \
              --chr-set {_MISS_BC_CHR_SET} no-xy no-mt \
              --set-missing-var-ids @:# \
              --make-bed \
              --threads {threads} \
              --out {params.base} \
              >> {log} 2>&1
        """

rule plink_missing_GTmiss_breed_choice:
    input:
        bed = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bed",
        bim = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bim",
        fam = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.fam"
    params:
        base = f"{_MISS_BC_DIR}/all_GTmiss"
    output:
        imiss = f"{_MISS_BC_DIR}/all_GTmiss.imiss",
        lmiss = f"{_MISS_BC_DIR}/all_GTmiss.lmiss"
    log:
        "logs/breed_choice_missingness/plink_missing_GTmiss_breed_choice.log"
    threads: 4
    resources:
        mem_mb  = 20000,
        runtime = 180
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {_MISS_BC_DIR}
        plink \
          --bed {input.bed} --bim {input.bim} --fam {input.fam} \
          --chr-set {_MISS_BC_CHR_SET} no-xy no-mt \
          --missing \
          --threads {threads} \
          --out {params.base} \
          >> {log} 2>&1
        """

rule plink_missing_reps_breed_choice:
    input:
        bed = f"{_MISS_BC_PLINK_DIR}/all_reps.bed",
        bim = f"{_MISS_BC_PLINK_DIR}/all_reps.bim",
        fam = f"{_MISS_BC_PLINK_DIR}/all_reps.fam"
    params:
        base = f"{_MISS_BC_DIR}/all_reps"
    output:
        imiss = f"{_MISS_BC_DIR}/all_reps.imiss",
        lmiss = f"{_MISS_BC_DIR}/all_reps.lmiss"
    log:
        "logs/breed_choice_missingness/plink_missing_reps_breed_choice.log"
    threads: 4
    resources:
        mem_mb  = 20000,
        runtime = 180
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {_MISS_BC_DIR}
        plink \
          --bed {input.bed} --bim {input.bim} --fam {input.fam} \
          --chr-set {_MISS_BC_CHR_SET} no-xy no-mt \
          --missing \
          --threads {threads} \
          --out {params.base} \
          >> {log} 2>&1
        """

rule plink_missing_matrix_thin_breed_choice:
    input:
        bed = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bed",
        bim = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.bim",
        fam = f"{_MISS_BC_PLINK_DIR}/all_GTmiss.fam"
    params:
        thin_base = f"{_MISS_BC_DIR}/all_GTmiss_thin50k",
        n_snps    = 50000
    output:
        traw = f"{_MISS_BC_DIR}/all_GTmiss_thin50k.traw"
    log:
        "logs/breed_choice_missingness/plink_missing_matrix_thin_breed_choice.log"
    threads: 2
    resources:
        mem_mb  = 12000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {_MISS_BC_DIR}
        plink \
          --bed {input.bed} --bim {input.bim} --fam {input.fam} \
          --chr-set {_MISS_BC_CHR_SET} no-xy no-mt \
          --thin-count {params.n_snps} \
          --recode A-transpose \
          --threads {threads} \
          --out {params.thin_base} \
          >> {log} 2>&1
        """

# Aggregate target: snakemake missingness_breed_choice
rule missingness_breed_choice:
    input:
        f"{_MISS_BC_DIR}/all_GTmiss.imiss",
        f"{_MISS_BC_DIR}/all_GTmiss.lmiss",
        f"{_MISS_BC_DIR}/all_reps.imiss",
        f"{_MISS_BC_DIR}/all_reps.lmiss",
        f"{_MISS_BC_DIR}/all_GTmiss_thin50k.traw"

rule bcf_stats_LD:
    input:
        vcf = "results/breed_choice/gvcf_LD/{autosome}_GTmiss.LDprune_bcf.gvcf"
    output:
        vcf = "results/breed_choice/genomics_stats/filt/{autosome}_LDprune_bcf.stats"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools stats {input.vcf} > {output.vcf}"
