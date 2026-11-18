# PCA/ADMIXTURE/KING keep list: drop F_MISS > 0.20, then repeatedly drop the sample with the most
# ngsRelate theta > 0.125 partners until no related pair is left.

_REL_DIR = f"{REL_ANALYSIS_DIR}/relatedness"

THETA_MAX = 0.125  # expected kinship of 2nd-degree relatives
REL_MIND = 0.20  # same as MIND in plink_clean_old.smk
MANUAL_EXCLUDE = {
    "TG1111D2629": "Mongolian; extreme heterozygosity outlier, excluded in the previous keep list",
}


# ngsRelate result used for the keep list (popgen.smk's ngsrelate_run can compute a new one)
rule stage_ngsrelate_breed_choice:
    input:
        "data/ngsrelate/ngsrelate.res"
    output:
        f"{_REL_DIR}/ngsrelate.res"
    shell:
        "cp {input} {output}"


rule vcf_sample_order_breed_choice:
    input:
        vcf = "results/breed_choice/gathered_vcf/all_GTmiss.gvcf"
    output:
        f"{_REL_DIR}/all_GTmiss.samples.txt"
    log:
        "results/breed_choice/logs/relatedness/vcf_sample_order.log"
    resources:
        mem_mb = 2000,
        runtime = 30
    conda:
        "../../envs/bcftools.yml"
    shell:
        "bcftools query -l {input.vcf} > {output} 2> {log}"


rule relatedness_exclude_breed_choice:
    input:
        res = f"{_REL_DIR}/ngsrelate.res",
        samples = f"{_REL_DIR}/all_GTmiss.samples.txt",
        imiss = "results/breed_choice/genomics_stats/missingness/all_GTmiss.imiss"
    output:
        keep = f"{_REL_DIR}/unrelated_samples_breed_choice.txt",
        drop = f"{_REL_DIR}/excluded_samples_breed_choice.csv",
        pairs = f"{_REL_DIR}/related_pairs_breed_choice.csv"
    params:
        theta_max = THETA_MAX,
        mind = REL_MIND
    log:
        "results/breed_choice/logs/relatedness/relatedness_exclude.log"
    resources:
        mem_mb = 4000,
        runtime = 30
    run:
        import pandas as pd
        from collections import defaultdict

        names = [l.strip() for l in open(input.samples) if l.strip()]
        res = pd.read_csv(input.res, sep="\t", usecols=["a", "b", "theta"])
        assert res[["a", "b"]].max().max() < len(names), "ngsRelate indices exceed VCF sample count"
        idx = dict(enumerate(names))  # ngsRelate a/b are 0-based VCF sample indices
        res["sample1"] = res["a"].map(idx)
        res["sample2"] = res["b"].map(idx)

        imiss = pd.read_csv(input.imiss, sep=r"\s+", dtype={"IID": str})
        fmiss = dict(zip(imiss["IID"], imiss["F_MISS"]))
        mind_fail = sorted(s for s in names if fmiss.get(s, 0) > params.mind)

        pairs = res[(res["theta"] > params.theta_max)
                    & ~res["sample1"].isin(mind_fail) & ~res["sample2"].isin(mind_fail)]
        adj = defaultdict(set)
        for s1, s2 in zip(pairs["sample1"], pairs["sample2"]):
            adj[s1].add(s2)
            adj[s2].add(s1)
        n_involved = len(adj)

        related_drop = []
        while any(adj.values()):
            v = max(adj, key=lambda n: (len(adj[n]), fmiss.get(n, 0), n))  # ties: higher F_MISS, then name
            related_drop.append(v)
            for nbr in adj[v]:
                adj[nbr].discard(v)
            adj[v] = set()

        reasons = {s: f"F_MISS > {params.mind}" for s in mind_fail}
        reasons.update({s: f"related (theta > {params.theta_max})" for s in related_drop})
        for s, why in MANUAL_EXCLUDE.items():
            reasons.setdefault(s, f"manual: {why}")
        keep = sorted(set(names) - set(reasons))

        os.makedirs(os.path.dirname(output.keep), exist_ok=True)
        with open(output.keep, "w") as fh:
            fh.writelines(s + "\n" for s in keep)
        pd.DataFrame(sorted(reasons.items()), columns=["sample_alias", "reason"]).to_csv(output.drop, index=False)
        pairs.assign(sample1_dropped=pairs["sample1"].isin(related_drop),
                     sample2_dropped=pairs["sample2"].isin(related_drop)) \
             .sort_values("theta", ascending=False) \
             [["sample1", "sample2", "theta", "sample1_dropped", "sample2_dropped"]] \
             .to_csv(output.pairs, index=False)

        with open(log[0], "w") as lf:
            lf.write(
                f"{len(names)} samples; {len(mind_fail)} fail F_MISS > {params.mind}; "
                f"{len(pairs)} pairs with theta > {params.theta_max} among the rest, involving "
                f"{n_involved} samples; {len(related_drop)} dropped by greedy removal; "
                f"{len(MANUAL_EXCLUDE)} manual exclusions; {len(keep)} kept.\n"
            )
