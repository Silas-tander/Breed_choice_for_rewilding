# Per-sample genotype counts (SnpEff VCFs) and the metadata table with heterozygosity, missingness, FROH.

META_BASE_BREED_CHOICE = "samples.csv"
META_DIR_BREED_CHOICE  = "results/breed_choice/metadata"


rule genotype_counts_breed_choice:
    input:
        vcf = "results/breed_choice/snpeff/{autosome}.ann.gvcf"
    params:
        chr = "{autosome}",
        outdir = "results/breed_choice/genotype_diagnostics/"
    output:
        zygosi = "results/breed_choice/genotype_diagnostics/{autosome}_genotype_counts.csv"
    log:
        "results/breed_choice/logs/genotype_counts/{autosome}.log"
    threads: 2
    resources:
        mem_mb = 8000,
        runtime = 240
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/genotype_counts.py {params.chr} {params.outdir} {input.vcf} &> {log}"


rule stage_metadata_breed_choice:
    input:
        f"metadata/{META_BASE_BREED_CHOICE}"
    output:
        f"{META_DIR_BREED_CHOICE}/{META_BASE_BREED_CHOICE}"
    shell:
        "cp {input} {output}"


rule update_metadata_breed_choice:
    input:
        meta = f"{META_DIR_BREED_CHOICE}/{META_BASE_BREED_CHOICE}",
        counts = expand(
            "results/breed_choice/genotype_diagnostics/{autosome}_genotype_counts.csv",
            autosome=config["autosomes_only"],
        ),
        fai = "data/ref_genome/EquCab3.0_chrom.fa.fai",
        roh = "results/breed_choice/ROH/bcftools/AF_all/chr_FROHints_bcf_cM_AF.csv",
        # missing_frac from GTmiss (repeats removed); SnpEff counts keep repeats -> *_wReps
        imiss = "results/breed_choice/genomics_stats/missingness/all_GTmiss.imiss"
    output:
        f"{META_DIR_BREED_CHOICE}/metadata_breed_choice_updated.csv"
    log:
        "results/breed_choice/logs/update_metadata_breed_choice.log"
    resources:
        mem_mb = 4000,
        runtime = 30
    conda:
        "../../envs/slurm_v2.yaml"
    shell:
        "python scripts/python/update_metadata.py --meta {input.meta} --counts {input.counts} "
        "--fai {input.fai} --roh {input.roh} --imiss {input.imiss} --counts-suffix wReps --out {output} &> {log}"
