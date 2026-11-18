# PCA/ADMIXTURE chain of plink_clean_old.smk + popgen.smk on alternative keep lists ({relset}),
# same file names under results/breed_choice/{relset}/.

# To add a relset: add its keep list here and its targets to rule all
RELSET_KEEP = {
    "unrelated_ngsKING0.0884": "results/breed_choice/unrelated_ngsKING0.0884/relatedness/unrelated_samples_breed_choice.txt",
}
_RS = "results/breed_choice/{relset}"
_RS_PLINK = f"{_RS}/plink/{PREFIX}"

wildcard_constraints:
    relset = "|".join(r.replace(".", r"\.") for r in RELSET_KEEP),


NGSKING_MAX = 0.0884  # KING 2nd/3rd-degree boundary (2^-3.5)

# Procedure of relatedness_exclude_breed_choice.smk on ngsRelate's KING-robust column
rule relatedness_exclude_ngsKING_breed_choice:
    input:
        res = f"{REL_ANALYSIS_DIR}/relatedness/ngsrelate.res",
        samples = f"{REL_ANALYSIS_DIR}/relatedness/all_GTmiss.samples.txt",
        imiss = "results/breed_choice/genomics_stats/missingness/all_GTmiss.imiss"
    output:
        keep = RELSET_KEEP["unrelated_ngsKING0.0884"],
        drop = "results/breed_choice/unrelated_ngsKING0.0884/relatedness/excluded_samples_breed_choice.csv",
        pairs = "results/breed_choice/unrelated_ngsKING0.0884/relatedness/related_pairs_breed_choice.csv"
    params:
        column = "KING",
        cutoff = NGSKING_MAX,
        mind = REL_MIND,
        manual = MANUAL_EXCLUDE
    log:
        "results/breed_choice/logs/relatedness/relatedness_exclude_ngsKING.log"
    resources:
        mem_mb = 4000,
        runtime = 30
    run:
        import pandas as pd
        from collections import defaultdict

        names = [l.strip() for l in open(input.samples) if l.strip()]
        res = pd.read_csv(input.res, sep="\t", usecols=["a", "b", params.column])
        assert res[["a", "b"]].max().max() < len(names), "ngsRelate indices exceed VCF sample count"
        idx = dict(enumerate(names))
        res["sample1"] = res["a"].map(idx)
        res["sample2"] = res["b"].map(idx)

        imiss = pd.read_csv(input.imiss, sep=r"\s+", dtype={"IID": str})
        fmiss = dict(zip(imiss["IID"], imiss["F_MISS"]))
        mind_fail = sorted(s for s in names if fmiss.get(s, 0) > params.mind)

        pairs = res[(res[params.column] > params.cutoff)
                    & ~res["sample1"].isin(mind_fail) & ~res["sample2"].isin(mind_fail)]
        adj = defaultdict(set)
        for s1, s2 in zip(pairs["sample1"], pairs["sample2"]):
            adj[s1].add(s2)
            adj[s2].add(s1)
        n_involved = len(adj)

        related_drop = []
        while any(adj.values()):
            v = max(adj, key=lambda n: (len(adj[n]), fmiss.get(n, 0), n))
            related_drop.append(v)
            for nbr in adj[v]:
                adj[nbr].discard(v)
            adj[v] = set()

        reasons = {s: f"F_MISS > {params.mind}" for s in mind_fail}
        reasons.update({s: f"related (ngsRelate {params.column} > {params.cutoff})" for s in related_drop})
        for s, why in params.manual.items():
            reasons.setdefault(s, f"manual: {why}")
        keep = sorted(set(names) - set(reasons))

        os.makedirs(os.path.dirname(output.keep), exist_ok=True)
        with open(output.keep, "w") as fh:
            fh.writelines(s + "\n" for s in keep)
        pd.DataFrame(sorted(reasons.items()), columns=["sample_alias", "reason"]).to_csv(output.drop, index=False)
        pairs.assign(sample1_dropped=pairs["sample1"].isin(related_drop),
                     sample2_dropped=pairs["sample2"].isin(related_drop)) \
             .sort_values(params.column, ascending=False) \
             [["sample1", "sample2", params.column, "sample1_dropped", "sample2_dropped"]] \
             .to_csv(output.pairs, index=False)

        with open(log[0], "w") as lf:
            lf.write(
                f"{len(names)} samples; {len(mind_fail)} fail F_MISS > {params.mind}; "
                f"{len(pairs)} pairs with ngsRelate {params.column} > {params.cutoff} among the rest, involving "
                f"{n_involved} samples; {len(related_drop)} dropped by greedy removal; "
                f"{len(params.manual)} manual exclusions; {len(keep)} kept.\n"
            )


rule relset_subset_rename_vcf:
    input:
        vcf = VCF_IN,
        samples = lambda wc: RELSET_KEEP[wc.relset],
        chr_dir = "results/breed_choice/chr"
    output:
        vcf = f"{_RS}/gathered_vcf/{PREFIX}.rename.gvcf.gz",
        csi = f"{_RS}/gathered_vcf/{PREFIX}.rename.gvcf.gz.csi"
    params:
        names = "results/breed_choice/chr/chr.names.txt"
    log:
        "results/breed_choice/logs/popstruct_relset/{relset}/subset_rename.log"
    threads: 2
    resources:
        mem_mb = 12000,
        runtime = 720 // 2
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        bcftools view -S {input.samples} --force-samples -Ou {input.vcf} 2> {log} \
          | bcftools annotate --rename-chrs {params.names} -Oz -o {output.vcf} 2>> {log}
        bcftools index --csi {output.vcf} 2>> {log}
        """


rule relset_plink_mind_nomiss_prune:
    input:
        vcf = f"{_RS}/gathered_vcf/{PREFIX}.rename.gvcf.gz"
    output:
        bed = f"{_RS_PLINK}.pruned.bed",
        bim = f"{_RS_PLINK}.pruned.bim",
        fam = f"{_RS_PLINK}.pruned.fam",
        mind_fam = f"{_RS_PLINK}.mind.fam",
        nomiss_bim = f"{_RS_PLINK}.mind.nomiss.bim",
        # read by filtration_counts_breed_choice
        base_bim = f"{_RS_PLINK}.bim",
        base_fam = f"{_RS_PLINK}.fam",
        mind_bim = f"{_RS_PLINK}.mind.bim",
        nomiss_fam = f"{_RS_PLINK}.mind.nomiss.fam"
    params:
        base = _RS_PLINK,
        window = WINDOW, step = STEP, r2 = R2, mind = MIND
    log:
        "results/breed_choice/logs/popstruct_relset/{relset}/plink_mind_nomiss_prune.log"
    threads: 2
    resources:
        mem_mb = 20000,
        runtime = 240
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        P="--chr-set {CHR_SET} no-xy no-mt --allow-extra-chr --threads {threads}"
        plink --vcf {input.vcf} --const-fid --set-missing-var-ids @:# --make-bed $P --out {params.base} > {log} 2>&1
        plink --bfile {params.base} --mind {params.mind} --make-bed $P --out {params.base}.mind >> {log} 2>&1
        plink --bfile {params.base}.mind --geno 0 --make-bed $P --out {params.base}.mind.nomiss >> {log} 2>&1
        test -s {output.nomiss_bim}
        plink --bfile {params.base}.mind.nomiss --indep-pairwise {params.window} {params.step} {params.r2} $P --out {params.base} >> {log} 2>&1
        plink --bfile {params.base}.mind.nomiss --extract {params.base}.prune.in --make-bed $P --out {params.base}.pruned >> {log} 2>&1
        echo "samples after --mind: $(wc -l < {output.mind_fam}); SNPs after --geno 0: $(wc -l < {output.nomiss_bim}); after pruning: $(wc -l < {output.bim})" >> {log}
        """


rule relset_pca:
    input:
        bed = f"{_RS_PLINK}.pruned.bed",
        bim = f"{_RS_PLINK}.pruned.bim",
        fam = f"{_RS_PLINK}.pruned.fam"
    output:
        eigenvec = f"{_RS_PLINK}.pruned.pca.eigenvec",
        eigenval = f"{_RS_PLINK}.pruned.pca.eigenval",
        eigenvec_maf = f"{_RS_PLINK}.pruned.pca.maf.eigenvec",
        eigenval_maf = f"{_RS_PLINK}.pruned.pca.maf.eigenval",
        maf_bed = f"{_RS_PLINK}.pruned.maf.bed",
        maf_bim = f"{_RS_PLINK}.pruned.maf.bim",
        maf_fam = f"{_RS_PLINK}.pruned.maf.fam"
    params:
        pruned = f"{_RS_PLINK}.pruned",
        npcs = 50,
        maf = 0.05
    log:
        "results/breed_choice/logs/popstruct_relset/{relset}/pca.log"
    threads: 4
    resources:
        mem_mb = 20000,
        runtime = 120
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        P="--chr-set {CHR_SET} no-xy no-mt --threads {threads}"
        plink --bfile {params.pruned} --pca {params.npcs} header tabs $P --out {params.pruned}.pca > {log} 2>&1
        plink --bfile {params.pruned} --maf {params.maf} --pca {params.npcs} header tabs $P --out {params.pruned}.pca.maf >> {log} 2>&1
        plink --bfile {params.pruned} --maf {params.maf} --make-bed $P --out {params.pruned}.maf >> {log} 2>&1
        """


# {adm}: '' = LD-pruned set, '_maf' = its MAF >= 0.05 subset; runs in its own dir so relsets don't collide
rule relset_admixture:
    input:
        bed = lambda wc: f"results/breed_choice/{wc.relset}/plink/{PREFIX}.pruned{'.maf' if wc.adm else ''}.bed",
        bim = lambda wc: f"results/breed_choice/{wc.relset}/plink/{PREFIX}.pruned{'.maf' if wc.adm else ''}.bim",
        fam = lambda wc: f"results/breed_choice/{wc.relset}/plink/{PREFIX}.pruned{'.maf' if wc.adm else ''}.fam"
    output:
        f"{_RS}/admixture_results{{adm}}/{{n}}/log.out"
    wildcard_constraints:
        adm = r"(|_maf)",
        n = r"\d+"
    threads: 4
    resources:
        mem_mb = 8000,
        runtime = 2*720
    conda:
        "../../envs/admixture.yaml"
    shell:
        r"""
        set -euo pipefail
        out=$(dirname {output}); bed=$(realpath {input.bed})
        mkdir -p $out && cd $out
        admixture --cv=10 -B1000 $bed {wildcards.n} -j{threads} > log.out
        """


# All eigenvalues, for % of total variance (the --pca 50 runs keep only 50)
rule pca_all_eigenval_breed_choice:
    input:
        bed = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned.bed",
        bim = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned.bim",
        fam = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned.fam"
    output:
        eigenval = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned.pca.all.eigenval"
    wildcard_constraints:
        pcaset = "|".join(r.replace(".", r"\.") for r in ["unrelated_theta0.125", *RELSET_KEEP])
    params:
        base = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned",
        tmp = "results/breed_choice/{pcaset}/plink/all_bcf_exclude2nd_dom2.pruned.pca.all"
    log:
        "results/breed_choice/logs/popstruct_relset/{pcaset}/pca_all_eigenval.log"
    threads: 4
    resources:
        mem_mb = 20000,
        runtime = 60
    conda:
        "../../envs/plink.yml"
    shell:
        r"""
        set -euo pipefail
        n=$(wc -l < {input.fam})
        plink --bfile {params.base} --chr-set {CHR_SET} no-xy no-mt --pca $n --threads {threads} --out {params.tmp} > {log} 2>&1
        rm -f {params.tmp}.eigenvec
        """
