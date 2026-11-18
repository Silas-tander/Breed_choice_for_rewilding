# Uncallable bp per autosome (assembly N gaps ∪ GenMap mask): denominator of per-sample heterozygosity.

rule autosome_uncallable_bp_breed_choice:
    input:
        fa = "data/ref_genome/EquCab3.0_chrom.fa",
        fai = "data/ref_genome/EquCab3.0_chrom.fa.fai",
        mask = "data/ref_genome/genmap_mask.merge.bed"  # same mask as remove_genmap
    output:
        "results/breed_choice/genomics_stats/autosome_uncallable_bp.tsv"
    params:
        x_contig = "NC_009175.3"  # chrX
    log:
        "results/breed_choice/logs/autosome_uncallable_bp.log"
    resources:
        mem_mb = 8000,
        runtime = 60
    conda:
        "../../envs/allel.yml"
    shell:
        "python scripts/python/autosome_uncallable_bp.py {input.fa} {input.fai} {input.mask} "
        "{params.x_contig} {output} > {log} 2>&1"
