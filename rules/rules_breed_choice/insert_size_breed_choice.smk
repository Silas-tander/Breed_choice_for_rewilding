# Mean insert size (samtools stats) of the in-house mols*/modoc* samples.

INSERT_SIZE_SAMPLES = [
    f"{prefix}{n}"
    for prefix in ("mols", "modoc")
    for n in range(1, 6)
]

INSERT_SIZE_DIR = "results/breed_choice/genomics_stats/insert_size"

rule insert_size_summary_breed_choice:
    input:
        expand(f"{INSERT_SIZE_DIR}/{{sample}}.insert_size.txt", sample=INSERT_SIZE_SAMPLES)
    output:
        f"{INSERT_SIZE_DIR}/insert_size_summary.csv"
    run:
        import csv

        rows = []
        for sample, path in zip(INSERT_SIZE_SAMPLES, input):
            with open(path) as fh:
                line = fh.read().strip()
            # samtools stats prints no insert size for unpaired data
            value = line.split("\t")[-1].strip() if line else "NA"
            rows.append((sample, value))

        with open(output[0], "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample", "insert_size_average"])
            writer.writerows(sorted(rows))

rule bam_insert_size_breed_choice:
    input:
        bam = "data/bam/{sample}.bam"
    output:
        f"{INSERT_SIZE_DIR}/{{sample}}.insert_size.txt"
    wildcard_constraints:
        sample = "|".join(INSERT_SIZE_SAMPLES)
    log:
        "logs/breed_choice_insert_size/{sample}.log"
    threads: 4
    resources:
        mem_mb = 4000,
        runtime = 120
    conda:
        "../../envs/samtools.yml"
    shell:
        r"""
        exec > {log} 2>&1
        set -euo pipefail
        samtools stats -@ {threads} {input.bam} | grep "insert size average" > {output} || true
        """
