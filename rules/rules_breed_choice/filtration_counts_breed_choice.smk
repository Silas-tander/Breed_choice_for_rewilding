# Sites left after each filtration step: bcftools stats (VCF steps), .bim/.fam (PLINK steps).

_FC_STATS = "results/breed_choice/genomics_stats/filt"
_FC_AUTOS = [f"chr{i}" for i in range(1, 32)]

# (step, stats suffix, description); the GenMap VCF is deleted, so it is counted with the SNP step
_FC_VCF_STEPS = [
    ("raw",         "_raw.stats",          "Raw joint-called gVCF (GenotypeGVCFs), all samples"),
    ("site",        "_bcf_harsh.stats",    "Site filter: 15×N < INFO/DP < 30×N, INFO/AC > 0, ExcessHet ≥ 1e-6"),
    ("genmap_snps", "_final.stats",        "GenMap mappability mask + biallelic SNPs only"),
    ("repeats",     ".reps.stats",         "Repeat mask (repeats_Equcab3.bed)"),
    ("gtmiss",      "_GTmiss.stats",       "Genotype filters → missing (FMT/DP < 8 or > 90, allele balance < 0.05, GQ < 8)"),
    ("gtmiss_wreps", "_GTmiss_wReps.stats", "Side branch: genotype filters with repeats kept (FMT/DP < 8 or > 90; heterozygosity/ROH)"),
]

def _fc_plink_steps(d, label):
    b = f"{d}/{PREFIX}"
    return [
        (f"{label}:subset", f"{b}.bim",             f"{b}.fam",             "Subset to keep list (PLINK bed)"),
        (f"{label}:mind",   f"{b}.mind.bim",        f"{b}.mind.fam",        "--mind 0.20"),
        (f"{label}:geno",   f"{b}.mind.nomiss.bim", f"{b}.mind.nomiss.fam", "--geno 0"),
        (f"{label}:prune",  f"{b}.pruned.bim",      f"{b}.pruned.fam",      "--indep-pairwise 50 5 0.2 (PCA input)"),
        (f"{label}:maf",    f"{b}.pruned.maf.bim",  f"{b}.pruned.maf.fam",  "--maf 0.05 (PCA .maf / ADMIXTURE _maf input)"),
    ]

_FC_PLINK_STEPS = (
    _fc_plink_steps(f"{REL_ANALYSIS_DIR}/plink", "theta0.125")
    + _fc_plink_steps("results/breed_choice/unrelated_ngsKING0.0884/plink", "ngsKING0.0884")
)
for _g in GENO_VALUES:
    _b = f"{REL_ANALYSIS_DIR}/plink_geno{_g}/{PREFIX}"
    _FC_PLINK_STEPS += [
        (f"theta0.125_geno{_g}:geno",  f"{_b}.mind.geno.bim", None, f"--geno {_g} (instead of --geno 0)"),
        (f"theta0.125_geno{_g}:prune", f"{_b}.pruned.bim",    None, "--indep-pairwise 50 5 0.2 (PCA input)"),
    ]

rule filtration_counts_breed_choice:
    input:
        stats = [f"{_FC_STATS}/{c}{suf}" for _, suf, _ in _FC_VCF_STEPS for c in _FC_AUTOS],
        bims = [bim for _, bim, _, _ in _FC_PLINK_STEPS],
        fams = [fam for _, _, fam, _ in _FC_PLINK_STEPS if fam]
    output:
        tsv = "results/breed_choice/genomics_stats/filtration_counts.tsv",
        per_chr = "results/breed_choice/genomics_stats/filtration_counts_per_chr.tsv"
    params:
        vcf_steps = _FC_VCF_STEPS,
        plink_steps = _FC_PLINK_STEPS,
        stats_dir = _FC_STATS,
        autos = _FC_AUTOS
    resources:
        mem_mb = 4000,
        runtime = 30
    run:
        import pandas as pd

        def sn(path):
            out = {}
            for line in open(path):
                if line.startswith("SN\t"):
                    _, _, key, val = line.rstrip("\n").split("\t")
                    out[key.rstrip(":")] = int(val)
                elif line.startswith("# AF"):
                    break
            return out

        def nlines(path):
            with open(path, "rb") as fh:
                return sum(1 for _ in fh)

        per_chr, rows = [], []
        for step, suf, desc in params.vcf_steps:
            for c in params.autos:
                s = sn(f"{params.stats_dir}/{c}{suf}")
                per_chr.append(dict(step=step, chrom=c, records=s["number of records"],
                                    snps=s["number of SNPs"], indels=s["number of indels"],
                                    multiallelic_snp_sites=s["number of multiallelic SNP sites"],
                                    samples=s["number of samples"]))
            sub = [r for r in per_chr if r["step"] == step]
            rows.append(dict(stage="VCF (chr1–31)", step=step, description=desc,
                             snps=sum(r["snps"] for r in sub), indels=sum(r["indels"] for r in sub),
                             records=sum(r["records"] for r in sub),
                             multiallelic_snp_sites=sum(r["multiallelic_snp_sites"] for r in sub),
                             samples=sub[0]["samples"]))
        for step, bim, fam, desc in params.plink_steps:
            rows.append(dict(stage=f"PLINK ({step.split(':')[0]})", step=step, description=desc,
                             snps=nlines(bim), samples=nlines(fam) if fam else None))

        pd.DataFrame(per_chr).to_csv(output.per_chr, sep="\t", index=False)
        pd.DataFrame(rows).to_csv(output.tsv, sep="\t", index=False)
