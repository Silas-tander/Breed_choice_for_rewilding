# KING (--ibdseg, --kinship) on all samples passing --mind 0.20; unpruned, as KING expects.

KING_ALL_DIR = "results/breed_choice/relatedness_king/king_mind0.20"
_KING_ALL_BED = "results/breed_choice/gathered_vcf/plink_missing_breed_choice/all_GTmiss"
_KING_MIND_BED = f"{KING_ALL_DIR}/all_GTmiss.mind"


rule king_input_mind_breed_choice:
    input:
        bed = f"{_KING_ALL_BED}.bed",
        bim = f"{_KING_ALL_BED}.bim",
        fam = f"{_KING_ALL_BED}.fam"
    output:
        bed = f"{_KING_MIND_BED}.bed",
        bim = f"{_KING_MIND_BED}.bim",
        fam = f"{_KING_MIND_BED}.fam",
        irem = f"{_KING_MIND_BED}.irem"
    params:
        inbase = _KING_ALL_BED,
        outbase = _KING_MIND_BED,
        mind = REL_MIND  # relatedness_exclude_breed_choice.smk
    log:
        "results/breed_choice/logs/relatedness/king_input_mind.log"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 60
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p {KING_ALL_DIR}
        plink --bfile {params.inbase} --chr-set 31 no-xy no-mt --allow-extra-chr \
              --mind {params.mind} --make-bed --threads {threads} --out {params.outbase} > {log} 2>&1
        touch {output.irem}   # plink only writes .irem when it removes someone
        """


# --sexchr 32: with KING's human coding horse chr23-31 would be dropped
rule king_ibdseg_mind_breed_choice:
    input:
        bed = f"{_KING_MIND_BED}.bed",
        bim = f"{_KING_MIND_BED}.bim",
        fam = f"{_KING_MIND_BED}.fam"
    output:
        seg = f"{KING_ALL_DIR}/king_ibs.seg",
        segments = f"{KING_ALL_DIR}/king_ibs.segments.gz"
    params:
        prefix = f"{KING_ALL_DIR}/king_ibs"
    log:
        "results/breed_choice/logs/relatedness/king_ibdseg_mind.log"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 120
    conda:
        "../../envs/king.yml"
    shell:
        r"""
        set -euo pipefail
        king -b {input.bed} --ibdseg --sexchr 32 --cpus {threads} --prefix {params.prefix} > {log} 2>&1
        """


# FID = IID, so --kinship writes a between-family .kin0
rule king_kinship_mind_breed_choice:
    input:
        bed = f"{_KING_MIND_BED}.bed",
        bim = f"{_KING_MIND_BED}.bim",
        fam = f"{_KING_MIND_BED}.fam"
    output:
        kin0 = f"{KING_ALL_DIR}/king_kinship.kin0"
    params:
        prefix = f"{KING_ALL_DIR}/king_kinship"
    log:
        "results/breed_choice/logs/relatedness/king_kinship_mind.log"
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 120
    conda:
        "../../envs/king.yml"
    shell:
        r"""
        set -euo pipefail
        king -b {input.bed} --kinship --sexchr 32 --cpus {threads} --prefix {params.prefix} > {log} 2>&1
        """
