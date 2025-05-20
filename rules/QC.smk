# QC rule gathered!
import os

rule fastqc:
    input:
        R1 = "data/pe/{accession}_1.fastq.gz",
        R2 = "data/pe/{accession}_2.fastq.gz"
    output:
        "results/QC/fastqc/{accession}_1_fastqc.zip",
        "results/QC/fastqc/{accession}_1_fastqc.html",
        "results/QC/fastqc/{accession}_2_fastqc.zip",
        "results/QC/fastqc/{accession}_2_fastqc.html"
    params:
        outdir = "results/QC/fastqc"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30
    conda:
        "../envs/fastQC.yml"
    shell:
        "fastqc {input.R1} {input.R2} -o {params.outdir}"

rule nested_fastqc:
    input:
        R1=lambda wildcards: validate_file_exists(
            f"data/mustang_fq/{wildcards.group}/{wildcards.accession}_R1_001.fastq.gz", wildcards
        ),
        R2=lambda wildcards: validate_file_exists(
            f"data/mustang_fq/{wildcards.group}/{wildcards.accession}_R2_001.fastq.gz", wildcards
        )
    output:
        R1_zip="results/QC/fastqc/{group}/{accession}_R1_001_fastqc.zip",
        R1_html="results/QC/fastqc/{group}/{accession}_R1_001_fastqc.html",
        R2_zip="results/QC/fastqc/{group}/{accession}_R2_001_fastqc.zip",
        R2_html="results/QC/fastqc/{group}/{accession}_R2_001_fastqc.html"
    params:
        outdir=lambda wildcards: f"results/QC/fastqc/{wildcards.group}"
    threads: 1
    resources:
        mem_mb=8000,
        runtime=620
    conda:
        "../envs/fastQC.yml"
    shell:
        "fastqc {input.R1} {input.R2} -o {params.outdir}"

# rule multiqc_read1:
#     input:
#         R1 = expand("results/QC/fastqc/{accession}_1_fastqc.zip", accession = config["accession"]),
#         # R2 = expand("results/QC/fastqc/{accession}_2_fastqc.zip", accession = config["accession"]),
#         report1 = expand("results/QC/fastqc/{accession}_1_fastqc.html", accession = config["accession"]),
#         # report2 = expand("results/QC/fastqc/{accession}_2_fastqc.html", accession = config["accession"]),
#     output:
#         qc = 'results/QC/multiqc_output/multiqc_report1.html'
#     threads: 4
#     resources:
#         mem_mb = 8000,
#         runtime = 30
#     conda:
#         "../envs/multiQC.yml"
#     shell:
#         """
#         multiqc -f -o {output.qc} {input.report1}
#         """

# rule multiqc_read2:
#     input:
#         # R1 = expand("results/QC/fastqc/{accession}_1_fastqc.gz", accession = config["accession"]),
#         R2 = expand("results/QC/fastqc/{accession}_2_fastqc.gz", accession = config["accession"]),
#         # report1 = expand("results/QC/fastqc/{accession}_1_fastqc.html", accession = config["accession"]),
#         report2 = expand("results/QC/fastqc/{accession}_2_fastqc.html", accession = config["accession"]),
#     output:
#         qc = 'results/QC/multiqc_output/multiqc_report2.html',
#         outdiro = 'results/QC/multiqc_output/'
#     threads: 4
#     resources:
#         mem_mb = 8000,
#         runtime = 30
#     conda:
#         "../envs/multiQC.yml"
#     shell:
#         """
#         multiqc -f -o {output.outdiro} {input.report2}
#         """