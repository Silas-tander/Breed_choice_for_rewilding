# Binned mean INFO/DP of the raw gVCFs and binned SNP density of the filtered VCFs (config: depth_bin_window).

rule bin_depth_breed_choice:
    input:
        vcf = "data/gvcf/{autosome}.gvcf.gz"
    params:
        window = config.get("depth_bin_window", 100000),
        chrom  = "{autosome}"  # raw CHROM is a RefSeq accession
    output:
        tsv = "results/breed_choice/depth_binned/{autosome}.depth_binned.tsv"
    log:
        "results/breed_choice/logs/depth_binned_breed_choice/bin_depth_breed_choice/{autosome}.log"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 180
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        (
        echo -e "CHROM\tBIN_START\tBIN_END\tMEAN_DP\tN_SITES" > {output.tsv}
        bcftools query -f '%POS\t%INFO/DP\n' {input.vcf} | \
        awk -v OFS='\t' -v w={params.window} -v chrom={params.chrom} '
            {{
                bin = int(($1-1)/w)
                sum[bin]+=$2; n[bin]++
            }}
            END {{
                for (b in sum) print chrom, b*w+1, (b+1)*w, sum[b]/n[b], n[b]
            }}' | sort -k2,2n >> {output.tsv}
        ) 2> {log}
        """

rule gather_depth_binned_breed_choice:
    input:
        tsvs = lambda wildcards: expand(
            "results/breed_choice/depth_binned/{autosome}.depth_binned.tsv",
            autosome=config["autosomes_only"]
        )
    output:
        "results/breed_choice/depth_binned/all_autosomes.depth_binned.tsv"
    log:
        "results/breed_choice/logs/depth_binned_breed_choice/gather_depth_binned_breed_choice.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        (
        head -n1 {input.tsvs[0]} > {output}
        for f in {input.tsvs}; do
            tail -n +2 "$f" >> {output}
        done
        ) 2> {log}
        """


# Bins without SNPs are not written (= 0)
rule bin_snp_density_breed_choice:
    input:
        vcf = "results/breed_choice/gvcf_filt/{autosome}.GTmiss.gvcf"
    params:
        window = config.get("depth_bin_window", 100000),
        chrom  = "{autosome}"
    output:
        tsv = "results/breed_choice/depth_binned/{autosome}.snp_density_binned.tsv"
    log:
        "results/breed_choice/logs/depth_binned_breed_choice/bin_snp_density_breed_choice/{autosome}.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        (
        echo -e "CHROM\tBIN_START\tBIN_END\tN_SNPS" > {output.tsv}
        bcftools query -i 'TYPE="snp"' -f '%POS\n' {input.vcf} | \
        awk -v OFS='\t' -v w={params.window} -v chrom={params.chrom} '
            {{ n[int(($1-1)/w)]++ }}
            END {{
                for (b in n) print chrom, b*w+1, (b+1)*w, n[b]
            }}' | sort -k2,2n >> {output.tsv}
        ) 2> {log}
        """

rule gather_snp_density_binned_breed_choice:
    input:
        tsvs = lambda wildcards: expand(
            "results/breed_choice/depth_binned/{autosome}.snp_density_binned.tsv",
            autosome=config["autosomes_only"]
        )
    output:
        "results/breed_choice/depth_binned/all_autosomes.snp_density_binned.tsv"
    log:
        "results/breed_choice/logs/depth_binned_breed_choice/gather_snp_density_binned_breed_choice.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 60
    conda:
        "../../envs/bcftools.yml"
    shell:
        r"""
        set -euo pipefail
        (
        head -n1 {input.tsvs[0]} > {output}
        for f in {input.tsvs}; do
            tail -n +2 "$f" >> {output}
        done
        ) 2> {log}
        """
