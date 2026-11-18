# PCA chain of plink_clean_old.smk with --geno {geno} instead of --geno 0 (GENO_VALUES: Snakefile).

_GENO_IN   = f"{PLINK_DIR}/{PREFIX}.mind"  # post --mind bed of plink_clean_old.smk
_GENO_BASE = f"{REL_ANALYSIS_DIR}/plink_geno{{geno}}/{PREFIX}"

rule geno_prune_pca_breed_choice:
    input:
        bed = f"{_GENO_IN}.bed",
        bim = f"{_GENO_IN}.bim",
        fam = f"{_GENO_IN}.fam"
    output:
        geno_bim = f"{_GENO_BASE}.mind.geno.bim",
        pruned_bim = f"{_GENO_BASE}.pruned.bim",
        eigenvec = f"{_GENO_BASE}.pruned.pca.eigenvec",
        eigenval = f"{_GENO_BASE}.pruned.pca.eigenval",
        eigenvec_maf = f"{_GENO_BASE}.pruned.pca.maf.eigenvec",
        eigenval_maf = f"{_GENO_BASE}.pruned.pca.maf.eigenval",
        eigenval_all = f"{_GENO_BASE}.pruned.pca.all.eigenval"
    wildcard_constraints:
        geno = "|".join(g.replace(".", r"\.") for g in GENO_VALUES)
    params:
        inp = _GENO_IN,
        base = _GENO_BASE,
        window = WINDOW, step = STEP, r2 = R2,
        npcs = 50,
        maf = 0.05
    log:
        "results/breed_choice/logs/popstruct_geno/geno{geno}.log"
    threads: 4
    resources:
        mem_mb = 40000,
        runtime = 360
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        P="--chr-set {CHR_SET} no-xy no-mt --threads {threads}"
        plink --bfile {params.inp} --geno {wildcards.geno} --make-bed $P --out {params.base}.mind.geno > {log} 2>&1
        plink --bfile {params.base}.mind.geno --indep-pairwise {params.window} {params.step} {params.r2} $P --out {params.base} >> {log} 2>&1
        plink --bfile {params.base}.mind.geno --extract {params.base}.prune.in --make-bed $P --out {params.base}.pruned >> {log} 2>&1
        plink --bfile {params.base}.pruned --pca {params.npcs} header tabs $P --out {params.base}.pruned.pca >> {log} 2>&1
        plink --bfile {params.base}.pruned --maf {params.maf} --pca {params.npcs} header tabs $P --out {params.base}.pruned.pca.maf >> {log} 2>&1
        n=$(wc -l < {params.base}.pruned.fam)
        plink --bfile {params.base}.pruned --pca $n $P --out {params.base}.pruned.pca.all >> {log} 2>&1
        rm -f {params.base}.pruned.pca.all.eigenvec
        echo "SNPs in: $(wc -l < {input.bim}); after --geno {wildcards.geno}: $(wc -l < {output.geno_bim}); after pruning: $(wc -l < {output.pruned_bim})" >> {log}
        """
