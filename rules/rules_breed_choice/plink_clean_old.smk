# Keep-list subset -> PLINK -> --mind 0.20 -> --geno 0 -> LD prune -> PCA, MAF subset, missingness, KING.

PREFIX = "all_bcf_exclude2nd_dom2"

VCF_IN   = "results/breed_choice/gathered_vcf/all_GTmiss.gvcf"
SAMPLES  = f"{REL_ANALYSIS_DIR}/relatedness/unrelated_samples_breed_choice.txt"  # ngsRelate keep list (relatedness_exclude_breed_choice.smk)

VCF      = f"{REL_ANALYSIS_DIR}/gathered_vcf/{PREFIX}.gvcf.gz"
VCF_RE   = f"{REL_ANALYSIS_DIR}/gathered_vcf/{PREFIX}.rename.gvcf.gz"

PLINK_DIR   = f"{REL_ANALYSIS_DIR}/plink"

BASE        = f"{PLINK_DIR}/{PREFIX}"

PEDMAP_BASE = "results/breed_choice/gathered_vcf/plink/all_GTmiss"

CHR_SET = 31
MIND    = 0.20  # per-sample missingness

KING_DIR  = f"{REL_ANALYSIS_DIR}/king"

RAW_BASE    = BASE
MIND_BASE   = f"{PLINK_DIR}/{PREFIX}.mind"
PRUNE_PREF  = f"{PLINK_DIR}/{PREFIX}"
PRUNED_BASE = f"{PLINK_DIR}/{PREFIX}.pruned"
KING_BASE   = f"{KING_DIR}/king_related"

WINDOW = 50  # --indep-pairwise window / step / r2
STEP   = 5
R2     = 0.2


# Subset the GTmiss VCF to the keep list (rule name is historical)
rule exclude_prze:
    input:
        vcf = VCF_IN,
        samples = SAMPLES
    output:
        vcf = VCF,
        csi = f"{VCF}.csi"
    threads: 2
    resources:
        mem_mb = 12000,
        runtime = 720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail

        # Write compressed VCF so downstream rules are consistent
        bcftools view -S {input.samples} --force-samples \
          -Oz -o {output.vcf} {input.vcf}

        bcftools index --csi {output.vcf}
        """


rule rename_chr_plink:
    input:
        vcf = VCF,
        csi = f"{VCF}.csi",
        chr_dir = "results/breed_choice/chr"  # ordering only: needs stage_chr_names_breed_choice
    params:
        names = "results/breed_choice/chr/chr.names.txt"
    output:
        vcf_re = VCF_RE,
        csi    = f"{VCF_RE}.csi"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 120
    conda:
        "../../envs/bcftools.yml"
    shell:
        r'''
        set -euo pipefail
        bcftools annotate --rename-chrs {params.names} {input.vcf} -Oz -o {output.vcf_re}
        bcftools index --csi {output.vcf_re}
        '''


rule PLINK_ped_map:
    input:
        vcf = VCF
    output:
        mapa = f"{PEDMAP_BASE}.map",
        ped  = f"{PEDMAP_BASE}.ped"
    params:
        base = PEDMAP_BASE,
        chr_set = CHR_SET
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        txt = f"{PEDMAP_BASE}.pedmap.log"
    shell:
        r"""
        set -euo pipefail
        plink --vcf {input.vcf} --recode --allow-extra-chr \
              --chr-set {params.chr_set} no-xy no-mt \
              --threads {threads} --out {params.base} \
              >> {log.txt} 2>&1
        """


rule PLINK_make_bed:
    input:
        vcf = VCF_RE,
        csi = f"{VCF_RE}.csi"
    output:
        bed = f"{BASE}.bed",
        bim = f"{BASE}.bim",
        fam = f"{BASE}.fam"
    params:
        base = BASE
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        txt = f"{BASE}.bed.log"
    shell:
        r"""
        set -euo pipefail
        # --set-missing-var-ids: the VCF has no variant IDs ("."), and without unique
        # IDs the later --extract {{prune}}.prune.in matches every "." variant, so LD
        # pruning silently kept all SNPs
        plink --vcf {input.vcf} --allow-extra-chr --const-fid \
              --chr-set {CHR_SET} no-xy no-mt \
              --set-missing-var-ids @:# \
              --make-bed \
              --threads {threads} --out {params.base} \
              >> {log.txt} 2>&1
        """


rule plink_filter_samples:
    input:
        bed = f"{RAW_BASE}.bed",
        bim = f"{RAW_BASE}.bim",
        fam = f"{RAW_BASE}.fam"
    output:
        bed = f"{MIND_BASE}.bed",
        bim = f"{MIND_BASE}.bim",
        fam = f"{MIND_BASE}.fam"
    params:
        base = MIND_BASE,
        mind = MIND
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        f"{MIND_BASE}.log"
    shell:
        r"""
        set -euo pipefail
        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --chr-set {CHR_SET} no-xy no-mt \
          --mind {params.mind} \
          --make-bed \
          --threads {threads} \
          --out {params.base} \
          >> {log} 2>&1
        """


# --geno 0: drop sites with any missing genotype before pruning
rule plink_drop_missing_sites_before_prune:
    input:
        bed = f"{MIND_BASE}.bed",
        bim = f"{MIND_BASE}.bim",
        fam = f"{MIND_BASE}.fam"
    output:
        bed = f"{PRUNE_PREF}.mind.nomiss.bed",
        bim = f"{PRUNE_PREF}.mind.nomiss.bim",
        fam = f"{PRUNE_PREF}.mind.nomiss.fam"
    params:
        outpref = f"{PRUNE_PREF}.mind.nomiss"
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNE_PREF}.mind.nomiss.log"
    shell:
        r"""
        set -euo pipefail
        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --chr-set {CHR_SET} no-xy no-mt \
          --geno 0 \
          --make-bed \
          --threads {threads} \
          --out {params.outpref} \
          >> {log} 2>&1

        nsnps=$(wc -l < {output.bim} || echo 0)
        echo "Remaining SNPs after --geno 0: $nsnps" >> {log}
        if [ "$nsnps" -eq 0 ]; then
          echo "ERROR: No SNPs left after removing missing-genotype sites." >> {log}
          exit 1
        fi
        """


rule plink_ld_prune_and_extract_nomiss:
    input:
        bed = f"{PRUNE_PREF}.mind.nomiss.bed",
        bim = f"{PRUNE_PREF}.mind.nomiss.bim",
        fam = f"{PRUNE_PREF}.mind.nomiss.fam"
    output:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    params:
        base   = PRUNE_PREF,
        pruned = PRUNED_BASE,
        window = WINDOW,
        step   = STEP,
        r2     = R2
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNE_PREF}.prune_and_extract_nomiss.log"
    shell:
        r"""
        set -euo pipefail

        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --chr-set {CHR_SET} no-xy no-mt \
          --allow-extra-chr \
          --indep-pairwise {params.window} {params.step} {params.r2} \
          --threads {threads} \
          --out {params.base} \
          >> {log} 2>&1

        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --chr-set {CHR_SET} no-xy no-mt \
          --allow-extra-chr \
          --extract {params.base}.prune.in \
          --make-bed \
          --threads {threads} \
          --out {params.pruned} \
          >> {log} 2>&1
        """


rule plink_pca_on_pruned:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        eigenvec = f"{PRUNE_PREF}.pruned.pca.eigenvec",
        eigenval = f"{PRUNE_PREF}.pruned.pca.eigenval"
    params:
        pruned = PRUNED_BASE,
        pcaout = f"{PRUNE_PREF}.pruned.pca",
        npcs   = 50
    threads: 4
    resources:
        mem_mb = 20000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNE_PREF}.pruned.pca.log"
    shell:
        r"""
        set -euo pipefail

        plink \
          --bfile {params.pruned} \
          --chr-set {CHR_SET} no-xy no-mt \
          --pca {params.npcs} header tabs \
          --threads {threads} \
          --out {params.pcaout} \
          >> {log} 2>&1

        test -s "{output.eigenvec}"
        test -s "{output.eigenval}"
        """


rule plink_pca_on_pruned_maf:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        eigenvec = f"{PRUNE_PREF}.pruned.pca.maf.eigenvec",
        eigenval = f"{PRUNE_PREF}.pruned.pca.maf.eigenval"
    params:
        pruned = PRUNED_BASE,
        pcaout = f"{PRUNE_PREF}.pruned.pca.maf",
        npcs   = 50,
        maf    = 0.05
    threads: 4
    resources:
        mem_mb = 20000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNE_PREF}.pruned.pca.maf.log"
    shell:
        r"""
        set -euo pipefail

        plink \
          --bfile {params.pruned} \
          --chr-set {CHR_SET} no-xy no-mt \
          --maf {params.maf} \
          --pca {params.npcs} header tabs \
          --threads {threads} \
          --out {params.pcaout} \
          >> {log} 2>&1

        test -s "{output.eigenvec}"
        test -s "{output.eigenval}"
        """


# MAF >= 0.05 subset of the pruned set, input of ADMIXTURE_MAF
rule plink_maf_filter_pruned:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        bed = f"{PRUNED_BASE}.maf.bed",
        bim = f"{PRUNED_BASE}.maf.bim",
        fam = f"{PRUNED_BASE}.maf.fam"
    params:
        pruned  = PRUNED_BASE,
        outpref = f"{PRUNED_BASE}.maf",
        maf     = 0.05
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 60
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNED_BASE}.maf.log"
    shell:
        r"""
        set -euo pipefail
        plink \
          --bfile {params.pruned} \
          --chr-set {CHR_SET} no-xy no-mt \
          --maf {params.maf} \
          --make-bed \
          --threads {threads} \
          --out {params.outpref} \
          >> {log} 2>&1
        """


rule plink_missing_report_on_pruned:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        imiss = f"{PRUNE_PREF}.pruned.imiss",
        lmiss = f"{PRUNE_PREF}.pruned.lmiss"
    params:
        outpref = f"{PRUNE_PREF}.pruned"
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 60
    conda:
        "../../envs/plink.yml"
    log:
        f"{PRUNE_PREF}.pruned.missing.log"
    shell:
        r"""
        set -euo pipefail
        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --chr-set {CHR_SET} no-xy no-mt \
          --missing \
          --threads {threads} \
          --out {params.outpref} \
          >> {log} 2>&1
        """


# One FID in the .fam, so KING --related writes only the within-family .kin
rule king_related:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        kin = f"{KING_BASE}.kin"
    params:
        prefix = KING_BASE
    threads: 1
    resources:
        mem_mb = 20000,
        runtime = 180
    conda:
        "../../envs/king.yml"
    log:
        f"{KING_BASE}.log"
    shell:
        r"""
        set -euo pipefail
        king -b {input.bed} --related --sexchr 32 --prefix {params.prefix} >> {log} 2>&1
        """

rule king_kinship:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        kin = f"{KING_DIR}/king_kinship.kin"
    params:
        prefix = f"{KING_DIR}/king_kinship"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/king.yml"
    log:
        f"{KING_DIR}/king_kinship.log"
    shell:
        r"""
        set -euo pipefail
        king -b {input.bed} --kinship --sexchr 32 --prefix {params.prefix} >> {log} 2>&1
        """

rule king_ibs:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        ibs = f"{KING_DIR}/king_ibs.seg"
    params:
        prefix = f"{KING_DIR}/king_ibs"
    threads: 4
    resources:
        mem_mb = 40000,
        runtime = 120
    conda:
        "../../envs/king.yml"
    log:
        f"{KING_DIR}/king_ibs.log"
    shell:
        r"""
        set -euo pipefail
        king -b {input.bed} --ibdseg --sexchr 32 --prefix {params.prefix} >> {log} 2>&1
        """


rule bed_to_vcf:
    input:
        bed = f"{PRUNED_BASE}.bed",
        bim = f"{PRUNED_BASE}.bim",
        fam = f"{PRUNED_BASE}.fam"
    output:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.pruned.vcf"
    params:
        outprefix = "results/breed_choice/gathered_vcf/all_GTmiss.pruned",
        chrset = CHR_SET
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    log:
        "results/breed_choice/gathered_vcf/all_GTmiss.pruned.bed_to_vcf.log"
    shell:
        r"""
        set -euo pipefail

        plink \
          --bed {input.bed} \
          --bim {input.bim} \
          --fam {input.fam} \
          --allow-extra-chr \
          --chr-set {params.chrset} no-xy no-mt \
          --recode vcf-iid \
          --threads {threads} \
          --out {params.outprefix} \
          >> {log} 2>&1
        """
