# Download rule
import os

# rule wget:
#     input:
#         srr = "{accession}"
#     output:


rule get_fastq:
    params:
        srr = "{accession}"
    output:
        "data/fq/{accession}.fastq.gz"
    threads:
        1
    resources:
        mem_mb = 8000,
        runtime = 4*720
    conda: "../envs/sratools.yml"
    shell:
        "fasterq-dump --concatenate-reads {params.srr} -o - | gzip > {output}"

rule get_fastq_pe_gz:
    output:
        # the wildcard name must be accession, pointing to an SRA number
        "data/pe/{accession}_1.fastq.gz",
        "data/pe/{accession}_2.fastq.gz",
    log:
        "logs/pe/{accession}.gz.log"
    params:
        extra="--skip-technical"
    threads: 1 
    resources:
        mem_mb = 8000,
        runtime = 4*720
    wrapper:
        "v3.5.0/bio/sra-tools/fasterq-dump"

rule sam_dump:
    params:
        bamsesh = "{bamsession}"
    output:
        "data/bams/{bamsession}.bam"
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 2*720
    conda: "../envs/sratools.yml"
    shell:
        "sam-dump {params.bamsesh} | samtools view -bS > {output}"


rule unpack:
    input:
        tarfile="data/mustang_tar/{tarfile}.tar"
    output:
        touch("data/mustang_tar/{tarfile}.done")
    params:
        extract_dir="data/mustang_fq/"  # Directory for extracted files
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 6*720
    shell:
        """
        mkdir -p {params.extract_dir}
        tar -xzvf {input.tarfile} -C {params.extract_dir}
        touch {output}
        """
